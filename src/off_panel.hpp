#ifndef VG_OFF_PANEL_HPP_INCLUDED
#define VG_OFF_PANEL_HPP_INCLUDED

/** \file
 * Off-panel detection, stage A: finding the bases that recur in the reads but that no graph walk
 * spells, while a window of reads is in memory. Each surviving candidate later becomes an edit
 * site (`snv_edit_site`), genotyped, phased and written like any other site.
 */

#include <cstdint>
#include <functional>
#include <ostream>
#include <string>
#include <vector>

#include <vg/vg.pb.h>

#include "handle.hpp"

namespace vg {

using namespace std;

/// The floors a candidate must pass to be kept, set low so that the call bar can be fitted from
/// the dump.
struct EditCountParams {
    /// Reads below this mapping quality are not counted.
    int min_mapq = 5;
    /// Bases below this quality are not counted, as the REF or the ALT. A read with no base
    /// qualities is counted whole.
    int min_base_quality = 20;
    /// The fewest fragments carrying the ALT.
    size_t min_fragments = 3;
    /// The smallest share of the fragments covering the base that carry the ALT.
    double min_fraction = 0.15;
};

/// A substitution the reads make at one base of one node: the node's base `ref` read as `alt`.
/// The offset and both bases are along the node's forward strand. Fragments are counted once
/// however many of their reads cover the base, by read name, which mates share; a fragment whose
/// reads disagree about the base is counted in `discordant_fragments` and nowhere else.
struct EditCandidate {
    nid_t node = 0;
    uint32_t offset = 0;
    char ref = 'N';
    char alt = 'N';
    uint32_t alt_fragments = 0;
    uint32_t ref_fragments = 0;
    /// Fragments with a third base.
    uint32_t other_fragments = 0;
    uint32_t discordant_fragments = 0;
    /// The mean mapping quality of the reads carrying the ALT.
    float mean_alt_mapq = 0.0f;

    /// The share of the fragments covering the base, discordant ones aside, that carry the ALT.
    double fraction() const;

    bool operator<(const EditCandidate& other) const {
        if (node != other.node) return node < other.node;
        if (offset != other.offset) return offset < other.offset;
        return alt < other.alt;
    }
};

/// The call bar: which candidates become edit sites, and which of those are staged once
/// genotyped. Applied after the floors of `EditCountParams`.
struct EditCallBar {
    /// The fewest fragments carrying the ALT.
    size_t min_alt_fragments = 5;
    /// The smallest ALT fraction.
    double min_fraction = 0.25;
    /// The lowest GQ at which a site called with the ALT is staged.
    double min_gq = 0.0;

    bool admits(const EditCandidate& candidate) const {
        return candidate.alt_fragments >= min_alt_fragments
               && candidate.fraction() >= min_fraction;
    }
};

/// What became of one candidate that passed the call bar: where it was placed, and its direct-pass
/// call.
struct EditOutcome {
    nid_t node = 0;
    uint32_t offset = 0;
    char alt = 'N';
    /// "top" or "nested" where the site was genotyped, or why it was not.
    string placed;
    string site;
    vector<int> genotype;
    double gq = 0.0;
    /// Whether it was staged, to be phased and written like any other site.
    bool staged = false;
};

/// The outcomes, filed from the threads that genotyped them without a lock.
class EditOutcomes {
public:
    explicit EditOutcomes(size_t threads) : by_thread(threads) {}
    void add(EditOutcome&& outcome);
    /// Every outcome, by node, offset and ALT.
    vector<EditOutcome> sorted() const;
private:
    vector<vector<EditOutcome>> by_thread;
};

/// Count the substitutions that the reads `for_each_read` visits make on the nodes `owns`
/// accepts, and return those that pass `params`, sorted. A substitution that a node next to its
/// own spells, as the other branch of a one-base bubble does, is a panel allele the mapper did
/// not take, and is left out.
vector<EditCandidate> count_snv_edits(
    const HandleGraph& graph,
    const function<void(const function<void(const Alignment&)>&)>& for_each_read,
    const function<bool(nid_t)>& owns, const EditCountParams& params);

/// Visit the candidates in `sorted`, as `count_snv_edits` sorts them, that lie on the nodes
/// `nodes`, which are sorted and free of duplicates.
void for_each_candidate_on(const vector<EditCandidate>& sorted, const vector<nid_t>& nodes,
                           const function<void(const EditCandidate&)>& iteratee);

}

#endif
