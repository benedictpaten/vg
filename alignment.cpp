#include "alignment.hpp"
#include "stream.hpp"

namespace vg {

int sam_for_each(string& filename, function<void(Alignment&)> lambda) {

    samFile *in = hts_open(filename.c_str(), "r");
    if (in == NULL) return 0;
    bam_hdr_t *hdr = sam_hdr_read(in);
    bam1_t *b = bam_init1();
    while (sam_read1(in, hdr, b) >= 0) {
        Alignment a = bam_to_alignment(b);
        lambda(a);
    }
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(in);
    return 1;

}

void write_alignments(std::ostream& out, vector<Alignment>& buf) {
    function<Alignment(uint64_t)> lambda =
        [&buf] (uint64_t n) {
        return buf[n];
    };
    stream::write(cout, buf.size(), lambda);
}

short quality_char_to_short(char c) {
    return static_cast<short>(c) - 33;
}

char quality_short_to_char(short i) {
    return static_cast<char>(i + 33);
}

void alignment_quality_short_to_char(Alignment& alignment) {
    const string& quality = alignment.quality();
    string squality;
    std::transform(quality.begin(), quality.end(), squality.begin(),
                   [](char c) { return (char)quality_short_to_char((int8_t)c); });
    alignment.set_quality(squality);
}

void alignment_quality_char_to_short(Alignment& alignment) {
    const string& quality = alignment.quality();
    string squality;
    std::transform(quality.begin(), quality.end(), squality.begin(),
                   [](char c) { return (int8_t)quality_char_to_short((char)c); });
    alignment.set_quality(squality);
}

Alignment bam_to_alignment(bam1_t* b) {

    Alignment alignment;

    uint32_t lqname = b->core.l_qname;
    string name; name.resize(lqname);
    int32_t lqseq = b->core.l_qseq;
    string sequence; sequence.resize(lqseq);
    string quality; quality.resize(lqseq);

    char* nameptr = bam_get_qname(b);
    memcpy((char*)name.c_str(), nameptr, lqname);

    uint8_t* qualptr = bam_get_qual(b);
    memcpy((char*)quality.c_str(), qualptr, lqseq);

    uint8_t* seqptr = bam_get_seq(b);
    for (int i = 0; i < lqseq; ++i) {
        sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
    }

    alignment.set_name(name);
    alignment.set_sequence(sequence);
    alignment.set_quality(quality);
    return alignment;
}

}