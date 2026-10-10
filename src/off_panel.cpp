#include "off_panel.hpp"

#include <algorithm>
#include <cctype>
#include <tuple>
#include <unordered_set>

#include <omp.h>

#include "utility.hpp"

namespace vg {

using namespace std;

double EditCandidate::fraction() const {
    const uint32_t covering = alt_fragments + ref_fragments + other_fragments;
    return covering == 0 ? 0.0 : (double)alt_fragments / (double)covering;
}

namespace {

/// One read base on a node base: `base` along the node's forward strand, or 0 where the read
/// matches the node.
struct Observation {
    nid_t node;
    uint32_t offset;
    char base;
    uint64_t fragment;
    int32_t mapq;

    bool operator<(const Observation& other) const {
        return tie(node, offset, base, fragment) < tie(other.node, other.offset, other.base,
                                                        other.fragment);
    }
};

/// Walk the read's aligned bases on the nodes `wanted` accepts, each of quality at least `min_q`,
/// giving each one's node, forward offset and forward base (0 for a match).
void for_each_aligned_base(const HandleGraph& graph, const Alignment& aln, int min_q,
                           const function<bool(nid_t)>& wanted,
                           const function<void(nid_t, uint32_t, char)>& visit) {
    const string& qual = aln.quality();
    size_t read_pos = 0;
    for (const Mapping& m : aln.path().mapping()) {
        const nid_t node = m.position().node_id();
        const bool take = wanted(node);
        const bool reverse = m.position().is_reverse();
        const size_t node_len = take ? graph.get_length(graph.get_handle(node)) : 0;
        size_t node_pos = (size_t)m.position().offset();
        for (const Edit& e : m.edit()) {
            const size_t from = (size_t)e.from_length();
            const size_t to = (size_t)e.to_length();
            if (take && from == to && from > 0) {
                for (size_t i = 0; i < from; ++i) {
                    const size_t r = read_pos + i;
                    if (!qual.empty() && (r >= qual.size() || (int)(uint8_t)qual[r] < min_q)) {
                        continue;
                    }
                    const size_t visited = node_pos + i;
                    if (visited >= node_len) {
                        continue;
                    }
                    const uint32_t forward = (uint32_t)(reverse ? node_len - 1 - visited : visited);
                    char base = 0;
                    if (!e.sequence().empty()) {
                        base = (char)toupper(e.sequence()[i]);
                        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
                            continue;
                        }
                        if (reverse) {
                            base = reverse_complement(base);
                        }
                    }
                    visit(node, forward, base);
                }
            }
            node_pos += from;
            read_pos += to;
        }
    }
}

/// Whether a node beside `node`, sharing a neighbour with it, spells the node with the base at
/// `offset` replaced by `alt`.
bool spelled_by_sibling(const HandleGraph& graph, nid_t node, uint32_t offset, char alt) {
    const handle_t handle = graph.get_handle(node);
    string wanted = graph.get_sequence(handle);
    wanted[offset] = alt;
    bool found = false;
    for (bool go_left : {true, false}) {
        graph.follow_edges(handle, go_left, [&](const handle_t& neighbour) {
            graph.follow_edges(neighbour, !go_left, [&](const handle_t& sibling) {
                if (graph.get_id(sibling) != node && graph.get_length(sibling) == wanted.size()
                    && graph.get_sequence(graph.get_is_reverse(sibling) ? graph.flip(sibling)
                                                                        : sibling) == wanted) {
                    found = true;
                }
                return !found;
            });
            return !found;
        });
    }
    return found;
}

}

vector<EditCandidate> count_snv_edits(
    const HandleGraph& graph,
    const function<void(const function<void(const Alignment&)>&)>& for_each_read,
    const function<bool(nid_t)>& owns, const EditCountParams& params) {

    vector<EditCandidate> kept;
    std::hash<string> name_hash;

    // The substitutions, by fragment.
    vector<Observation> alts;
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq) {
            return;
        }
        const uint64_t fragment = name_hash(aln.name());
        for_each_aligned_base(graph, aln, params.min_base_quality, owns,
                              [&](nid_t node, uint32_t offset, char base) {
                                  if (base != 0) {
                                      alts.push_back(Observation{node, offset, base, fragment,
                                                                 aln.mapping_quality()});
                                  }
                              });
    });
    if (alts.empty()) {
        return kept;
    }
    std::sort(alts.begin(), alts.end());

    // The (node, offset) of each substitution that enough fragments carry, before the fragments
    // covering the base are counted.
    vector<pair<nid_t, uint32_t>> positions;
    unordered_set<nid_t> position_nodes;
    for (size_t i = 0; i < alts.size();) {
        size_t j = i;
        size_t fragments = 0;
        for (; j < alts.size() && alts[j].node == alts[i].node && alts[j].offset == alts[i].offset
               && alts[j].base == alts[i].base;
             ++j) {
            if (j == i || alts[j].fragment != alts[j - 1].fragment) {
                ++fragments;
            }
        }
        if (fragments >= params.min_fragments) {
            positions.emplace_back(alts[i].node, alts[i].offset);
            position_nodes.insert(alts[i].node);
        }
        i = j;
    }
    if (positions.empty()) {
        return kept;
    }
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());

    // Every base the reads put at those positions, by fragment.
    vector<Observation> covering;
    for_each_read([&](const Alignment& aln) {
        if (aln.mapping_quality() < params.min_mapq) {
            return;
        }
        const uint64_t fragment = name_hash(aln.name());
        for_each_aligned_base(
            graph, aln, params.min_base_quality,
            [&](nid_t node) { return owns(node) && position_nodes.count(node); },
            [&](nid_t node, uint32_t offset, char base) {
                if (std::binary_search(positions.begin(), positions.end(),
                                       make_pair(node, offset))) {
                    covering.push_back(Observation{node, offset, base, fragment, 0});
                }
            });
    });
    std::sort(covering.begin(), covering.end(),
              [](const Observation& a, const Observation& b) {
                  return tie(a.node, a.offset, a.fragment, a.base)
                         < tie(b.node, b.offset, b.fragment, b.base);
              });

    for (size_t i = 0; i < covering.size();) {
        // One position.
        size_t end = i;
        while (end < covering.size() && covering[end].node == covering[i].node
               && covering[end].offset == covering[i].offset) {
            ++end;
        }
        // Each fragment's one base, or 1 for a fragment whose reads disagree.
        vector<char> fragment_base;
        for (size_t f = i; f < end;) {
            size_t g = f;
            bool agree = true;
            while (g < end && covering[g].fragment == covering[f].fragment) {
                agree = agree && covering[g].base == covering[f].base;
                ++g;
            }
            fragment_base.push_back(agree ? covering[f].base : (char)1);
            f = g;
        }
        const nid_t node = covering[i].node;
        const uint32_t offset = covering[i].offset;
        const char ref = (char)toupper(graph.get_base(graph.get_handle(node), offset));
        for (char alt : {'A', 'C', 'G', 'T'}) {
            if (alt == ref) {
                continue;
            }
            EditCandidate c;
            c.node = node;
            c.offset = offset;
            c.ref = ref;
            c.alt = alt;
            for (char b : fragment_base) {
                if (b == 1) {
                    ++c.discordant_fragments;
                } else if (b == 0) {
                    ++c.ref_fragments;
                } else if (b == alt) {
                    ++c.alt_fragments;
                } else {
                    ++c.other_fragments;
                }
            }
            if (c.alt_fragments < params.min_fragments || c.fraction() < params.min_fraction) {
                continue;
            }
            if (spelled_by_sibling(graph, node, offset, alt)) {
                continue;
            }
            // The mean mapping quality of the reads carrying the ALT.
            Observation probe{node, offset, alt, 0, 0};
            double total = 0.0;
            size_t n = 0;
            for (auto it = std::lower_bound(alts.begin(), alts.end(), probe);
                 it != alts.end() && it->node == node && it->offset == offset && it->base == alt;
                 ++it) {
                total += it->mapq;
                ++n;
            }
            c.mean_alt_mapq = n == 0 ? 0.0f : (float)(total / (double)n);
            kept.push_back(c);
        }
        i = end;
    }
    std::sort(kept.begin(), kept.end());
    return kept;
}

void EditOutcomes::add(EditOutcome&& outcome) {
    by_thread.at(omp_get_thread_num()).push_back(std::move(outcome));
}

vector<EditOutcome> EditOutcomes::sorted() const {
    vector<EditOutcome> all;
    for (const auto& list : by_thread) {
        all.insert(all.end(), list.begin(), list.end());
    }
    std::sort(all.begin(), all.end(), [](const EditOutcome& a, const EditOutcome& b) {
        return tie(a.node, a.offset, a.alt) < tie(b.node, b.offset, b.alt);
    });
    return all;
}

void for_each_candidate_on(const vector<EditCandidate>& sorted, const vector<nid_t>& nodes,
                           const function<void(const EditCandidate&)>& iteratee) {
    auto it = sorted.begin();
    for (nid_t node : nodes) {
        it = std::lower_bound(it, sorted.end(), node,
                              [](const EditCandidate& c, nid_t n) { return c.node < n; });
        for (; it != sorted.end() && it->node == node; ++it) {
            iteratee(*it);
        }
    }
}

}
