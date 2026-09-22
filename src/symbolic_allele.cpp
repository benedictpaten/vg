#include "symbolic_allele.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>

#include <unordered_map>

namespace vg {

/// Where each node id appears in the traversal, so "does this chain close later" is a lookup
/// rather than a rescan. Traversals of a large snarl run to thousands of visits and every visit
/// asks the question.
static unordered_map<nid_t, vector<int>> index_positions(const SnarlTraversal& trav) {
    unordered_map<nid_t, vector<int>> at;
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& v = trav.visit(i);
        if (!v.has_snarl()) {
            at[v.node_id()].push_back(i);
        }
    }
    return at;
}

/// The chain a child snarl belongs to, as its boundary node ids. Falls back to the snarl's own
/// boundaries when it is not in a chain the manager knows about, which is the trivial-chain case:
/// chains_of() presents unary snarls and trivial-chain snarls as their own chains anyway, so the
/// fallback agrees with the general case rather than being a special one.
static pair<nid_t, nid_t> chain_bounds(const Snarl* child, const SnarlManager& snarl_manager) {
    const Chain* chain = snarl_manager.chain_of(child);
    if (chain != nullptr && !chain->empty()) {
        return make_pair(get_start_of(*chain).node_id(), get_end_of(*chain).node_id());
    }
    return make_pair(child->start().node_id(), child->end().node_id());
}

/// The snarl `site` refers to, as the manager knows it, or null when the two disagree. Shared by
/// projection and by `symbolic_site_resolvable` so the two can never drift apart.
static const Snarl* resolve_site(const Snarl& site, const SnarlManager& snarl_manager,
                                 bool* out_reversed = nullptr) {
    if (out_reversed != nullptr) {
        *out_reversed = false;
    }
    const Snarl* site_ptr = snarl_manager.into_which_snarl(site.start().node_id(),
                                                          site.start().backward());
    if (site_ptr == nullptr) {
        return nullptr;
    }
    // The point of this test is to confirm we got back *this* snarl rather than some other snarl
    // reachable by entering that node, so it compares boundary nodes -- but it must accept them in
    // EITHER order, because a snarl and its reversal are the same snarl.
    //
    // `flip_snarl` (graph_caller.cpp) reverses a snarl whose reference path runs backwards, and the
    // caller then works on that reversed copy: its start node is the original END node. Requiring
    // start-to-start therefore failed for every such snarl, `site_ptr` came back null, `is_child`
    // was false at every visit, and the projection degenerated to a bare node list with no chain
    // symbols at all -- silently turning symbolic collapsing off for 7.4% of chr20's sites, a
    // feature worth SNV F1 0.9752 -> 0.9833 where it does run.
    //
    // Accepting the reversed pairing does not weaken the test. The boundary index maps a snarl's
    // start and its reversed end to the same snarl (SnarlManager::snarl_boundary_index), so a
    // reversed match identifies the same snarl by the same two nodes, just entered from the other
    // side. Everything downstream is orientation-independent: `parent_of` compares canonical
    // pointers, `chain_bounds` returns node ids, and a chain symbol's direction is taken from which
    // boundary the traversal meets first rather than from the site's own orientation.
    auto same_visit = [](const Visit& a, const Visit& b) {
        return a.node_id() == b.node_id() && a.backward() == b.backward();
    };
    // Node ids only, exactly as before this change: the forward branch is pre-existing behaviour and
    // tightening it would reject sites that resolve today, which is not this fix's business.
    const bool forward = site_ptr->start().node_id() == site.start().node_id() &&
                         site_ptr->end().node_id() == site.end().node_id();
    // Full visits, orientation included. This branch is the only acceptance this change adds, so it
    // is held to the exact test rather than inheriting the forward branch's laxity: a reversed snarl
    // is `reverse`-and-swap of the canonical one, and nothing less specific than that should pass.
    const bool reversed = same_visit(site_ptr->start(), reverse(site.end())) &&
                          same_visit(site_ptr->end(), reverse(site.start()));
    if (!forward && !reversed) {
        return nullptr;
    }
    if (out_reversed != nullptr) {
        *out_reversed = reversed && !forward;
    }
    return site_ptr;
}

bool symbolic_site_resolvable(const Snarl& site, const SnarlManager& snarl_manager,
                              bool* out_reversed) {
    return resolve_site(site, snarl_manager, out_reversed) != nullptr;
}

pair<nid_t, nid_t> chain_bounds_of(const Snarl* child, const SnarlManager& snarl_manager) {
    return chain_bounds(child, snarl_manager);
}

SymbolicAllele symbolic_allele(const SnarlTraversal& trav, const Snarl& site,
                               const SnarlManager& snarl_manager,
                               vector<pair<int, int>>* out_visit_ranges) {
    SymbolicAllele out;
    if (out_visit_ranges != nullptr) {
        out_visit_ranges->clear();
    }
    unordered_map<nid_t, vector<int>> at = index_positions(trav);

    // The snarl we are projecting, as the manager knows it, so a candidate child can be tested for
    // being a genuine child of *this* site rather than merely a snarl the traversal walks into.
    //
    // Null here means no child is ever recognised and the projection degenerates to the plain node
    // list. That is not hypothetical: `flip_snarl` reverses a snarl whose reference path runs
    // backwards, and the reversed copy's start node is the original END node, so the identity test
    // below fails and symbolic collapsing is inert for the whole snarl. See
    // `symbolic_site_resolvable`, which is how that population is counted.
    const Snarl* site_ptr = resolve_site(site, snarl_manager);

    int i = 0;
    while (i < trav.visit_size()) {
        const Visit& v = trav.visit(i);

        // A visit already carrying a Snarl is a symbol as it stands.
        if (v.has_snarl()) {
            SymbolicStep step;
            step.id = v.snarl().start().node_id();
            step.end_id = v.snarl().end().node_id();
            step.backward = v.backward();
            out.push_back(step);
            if (out_visit_ranges != nullptr) {
                out_visit_ranges->emplace_back(i, i + 1);
            }
            ++i;
            continue;
        }

        const nid_t node = v.node_id();
        // The snarl entered by traversing into this node in this orientation, if any.
        const Snarl* child = snarl_manager.into_which_snarl(node, v.backward());
        bool symbolised = false;

        // Only a genuine child of this site may be symbolised. Comparing the *chain's* boundaries
        // against the site's is not enough: a site that is itself a member of a longer chain sees
        // that enclosing chain's bounds, which differ from its own, and collapses its own interior
        // into one symbol -- making every allele equal and silently erasing the variant. A snarl
        // (5,9) inside a chain spanning 2..9 did exactly that.
        bool is_child = child != nullptr && site_ptr != nullptr &&
                        snarl_manager.parent_of(child) == site_ptr;
        if (is_child) {
            pair<nid_t, nid_t> bounds = chain_bounds(child, snarl_manager);
            {
                // Leave the chain at whichever of its boundaries this traversal reaches next; a
                // chain can be crossed in either direction, so both are candidates.
                int exit = -1;
                for (nid_t boundary : {bounds.second, bounds.first}) {
                    if (boundary == node) {
                        continue;
                    }
                    auto found = at.find(boundary);
                    if (found == at.end()) {
                        continue;
                    }
                    // index_positions fills each vector in increasing visit order, so it is sorted
                    // and the first entry past i is the nearest exit. Scanning the whole vector
                    // instead is O(k) in the node's recurrence, which is the dominant cost in
                    // satellite where one node recurs thousands of times in a single traversal.
                    auto it = std::upper_bound(found->second.begin(), found->second.end(), i);
                    if (it != found->second.end() && (exit < 0 || *it < exit)) {
                        exit = *it;
                    }
                }
                if (exit > i) {
                    SymbolicStep step;
                    step.id = bounds.first;
                    step.end_id = bounds.second;
                    // The chain is traversed backward when its recorded end is met before its start.
                    step.backward = (node == bounds.second);
                    out.push_back(step);
                    if (out_visit_ranges != nullptr) {
                        out_visit_ranges->emplace_back(i, exit);
                    }
                    // Resume *at* the exit boundary, which belongs to both the chain and whatever
                    // follows it, so it is not consumed.
                    i = exit;
                    symbolised = true;
                }
                // A chain entered and never left within this traversal falls through deliberately
                // and is emitted as a plain node. Swallowing the remainder would drop real
                // differences and make unrelated alleles compare equal.
            }
        }

        if (!symbolised) {
            SymbolicStep step;
            step.id = node;
            step.backward = v.backward();
            out.push_back(step);
            if (out_visit_ranges != nullptr) {
                out_visit_ranges->emplace_back(i, i + 1);
            }
            ++i;
        }
    }
    return out;
}

bool symbolically_equal(const SnarlTraversal& a, const SnarlTraversal& b, const Snarl& site,
                        const SnarlManager& snarl_manager) {
    return symbolic_allele(a, site, snarl_manager) == symbolic_allele(b, site, snarl_manager);
}


vector<DiffBlock> symbolic_diff(const SymbolicAllele& ref, const SymbolicAllele& alt,
                                bool* out_degraded, vector<int>* out_alt_before_ref) {
    if (out_degraded != nullptr) {
        *out_degraded = false;
    }

    const size_t m = ref.size();
    const size_t n = alt.size();

    // Filled for every exit path, so a caller never reads a stale or short vector.
    auto trivial_map = [&](size_t consumed_before_end) {
        if (out_alt_before_ref == nullptr) {
            return;
        }
        out_alt_before_ref->assign(m + 1, 0);
        (*out_alt_before_ref)[m] = (int)consumed_before_end;
    };

    if (m == 0 && n == 0) {
        trivial_map(0);
        return {};
    }
    if (m == 0 || n == 0) {
        // Wholly an insertion or wholly a deletion. No alignment to compute, and not a degradation.
        trivial_map(n);
        return {DiffBlock{0, (int)m, 0, (int)n}};
    }

    // Ukkonen's band, not the full rectangle. Every cell on an optimal path satisfies
    // |i - j| <= D, where D is the edit distance, so a band of half-width k >= D holds the whole
    // optimum and everything outside it is provably worse. k starts at the length difference (the
    // minimum any alignment must spend) and doubles until the corner value certifies itself by
    // coming out <= k, which is Ukkonen's termination test. Work is O((m + n) * D) rather than
    // O(m * n), and for two alleles of one site D is the actual variation between them, not their
    // length.
    //
    // This REPLACES a 4M-cell cap that degraded large pairs to one whole-allele block. Nothing
    // degrades now, so `out_degraded` is always false; it is kept because callers read it.
    //
    // Out-of-band cells read as INF, which is what makes the traceback below identical to the full
    // matrix's: a cell outside the band has true value > D >= the value of any cell on the optimal
    // path, so `at(i,j) == at(i-1,j) + 1` could never have held for it. Reading INF fails the same
    // test the true value would have failed, so the same move is chosen at every step.
    const uint32_t INF = std::numeric_limits<uint32_t>::max() / 4;
    const size_t max_k = std::max(m, n);
    size_t band_k = std::max<size_t>(1, m > n ? m - n : n - m);
    vector<uint32_t> cost;
    size_t stride = 0;

    // Read: INF outside the band. Write: only ever called in band.
    auto cell = [&](size_t i, size_t j) -> uint32_t {
        long long off = (long long)j - (long long)i + (long long)band_k;
        if (off < 0 || off >= (long long)stride) {
            return INF;
        }
        return cost[i * stride + (size_t)off];
    };
    auto put = [&](size_t i, size_t j, uint32_t v) {
        cost[i * stride + (size_t)((long long)j - (long long)i + (long long)band_k)] = v;
    };

    while (true) {
        stride = 2 * band_k + 1;
        cost.assign(stride * (m + 1), INF);
        // Row and column initialisations, clipped to the band.
        for (size_t i = 0; i <= m && i <= band_k; ++i) {
            put(i, 0, (uint32_t)i);
        }
        for (size_t j = 0; j <= n && j <= band_k; ++j) {
            put(0, j, (uint32_t)j);
        }
        for (size_t i = 1; i <= m; ++i) {
            const size_t jlo = (i > band_k) ? i - band_k : 1;
            const size_t jhi = std::min(n, i + band_k);
            for (size_t j = jlo == 0 ? 1 : jlo; j <= jhi; ++j) {
                uint32_t diag = cell(i - 1, j - 1);
                uint32_t del = cell(i - 1, j);
                uint32_t ins = cell(i, j - 1);
                if (diag < INF) {
                    diag += (ref[i - 1] == alt[j - 1] ? 0u : 1u);
                }
                if (del < INF) {
                    del += 1u;
                }
                if (ins < INF) {
                    ins += 1u;
                }
                put(i, j, std::min(diag, std::min(del, ins)));
            }
        }
        if (cell(m, n) <= (uint32_t)band_k || band_k >= max_k) {
            break;
        }
        band_k = std::min(max_k, band_k * 2);
    }
    auto at = [&](size_t i, size_t j) -> uint32_t { return cell(i, j); };

    // Traceback. The preference order here IS the tie-break contract in the header: diagonal first
    // (match before substitute, since a match is the zero-cost diagonal), then deletion, then
    // insertion. Walked backwards, so the ops come out reversed and are flipped below.
    enum Op { OP_MATCH, OP_SUB, OP_DEL, OP_INS };
    vector<Op> ops;
    ops.reserve(m + n);
    {
        size_t i = m;
        size_t j = n;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0) {
                bool equal = ref[i - 1] == alt[j - 1];
                if (at(i, j) == at(i - 1, j - 1) + (equal ? 0u : 1u)) {
                    ops.push_back(equal ? OP_MATCH : OP_SUB);
                    --i;
                    --j;
                    continue;
                }
            }
            if (i > 0 && at(i, j) == at(i - 1, j) + 1u) {
                ops.push_back(OP_DEL);
                --i;
                continue;
            }
            // j > 0 necessarily: the row and column initialisations make the remaining move legal.
            ops.push_back(OP_INS);
            --j;
        }
    }
    std::reverse(ops.begin(), ops.end());

    if (out_alt_before_ref != nullptr) {
        // Entry i is the alt index on arrival at reference step i, which is after everything that
        // consumed reference steps below i but before anything inserted at boundary i. That is the
        // convention the join needs: an insertion at a boundary belongs to the block that owns the
        // boundary, not to the reference step after it.
        out_alt_before_ref->assign(m + 1, 0);
        size_t ri = 0;
        size_t ai = 0;
        for (Op op : ops) {
            if (op == OP_INS) {
                ++ai;
            } else if (op == OP_DEL) {
                ++ri;
                (*out_alt_before_ref)[ri] = (int)ai;
            } else {
                ++ri;
                ++ai;
                (*out_alt_before_ref)[ri] = (int)ai;
            }
        }
    }

    // Aggregate every maximal run of non-match ops into one block. A substitution is a difference,
    // so it joins the run it sits in rather than splitting it -- which is what makes an interior
    // mismatch inside a longer difference come out as one record instead of three.
    vector<DiffBlock> out;
    size_t ri = 0;
    size_t ai = 0;
    size_t k = 0;
    while (k < ops.size()) {
        if (ops[k] == OP_MATCH) {
            ++ri;
            ++ai;
            ++k;
            continue;
        }
        DiffBlock block;
        block.ref_begin = (int)ri;
        block.alt_begin = (int)ai;
        while (k < ops.size() && ops[k] != OP_MATCH) {
            if (ops[k] == OP_SUB) {
                ++ri;
                ++ai;
            } else if (ops[k] == OP_DEL) {
                ++ri;
            } else {
                ++ai;
            }
            ++k;
        }
        block.ref_end = (int)ri;
        block.alt_end = (int)ai;
        out.push_back(block);
    }
    return out;
}

ostream& operator<<(ostream& out, const DiffBlock& block) {
    return out << "[ref " << block.ref_begin << "," << block.ref_end
               << " alt " << block.alt_begin << "," << block.alt_end << "]";
}

ostream& operator<<(ostream& out, const SymbolicAllele& allele) {
    for (const SymbolicStep& s : allele) {
        out << (s.backward ? '<' : '>');
        if (s.is_chain()) {
            out << "C" << s.id << "_" << s.end_id;
        } else {
            out << s.id;
        }
    }
    return out;
}

}
