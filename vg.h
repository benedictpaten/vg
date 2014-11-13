#ifndef VG_H
#define VG_H

#include <vector>
#include <set>
#include <string>
#include "gssw.h"
#include "vg.pb.h"
#include "Variant.h"
#include "Fasta.h"


namespace vg {


class VariantGraph : public Graph {

public:

    // nodes by id
    map<int64_t, Node*> node_by_id;

    // nodes by position in nodes repeated field
    // this is critical to allow fast deletion of nodes
    map<Node*, int> node_index;

    // edges indexed by nodes they connect
    map<int64_t, map<int64_t, Edge*> > edge_from_to;
    map<int64_t, map<int64_t, Edge*> > edge_to_from;

    // edges by position in edges repeated field
    // same as for nodes, this allows fast deletion
    map<Edge*, int> edge_index;

    // constructors
    //VariantGraph(void) { };
    // construct from protobufs
    VariantGraph(istream& in);
    VariantGraph(Graph& graph);
    VariantGraph(vector<Node>& nodes);

    // construct from VCF records
    VariantGraph(vcf::VariantCallFile& variantCallFile, FastaReference& reference);

    // use the VariantGraph class to generate ids
    Node* create_node(string seq);
    void destroy_node(Node* node);

    Edge* create_edge(Node* from, Node* to);
    Edge* create_edge(int64_t from, int64_t to);
    void destroy_edge(Edge* edge);

    // utilities
    void divide_node(Node* node, int pos, Node*& left, Node*& right);
    void divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right);
    //void node_replace_prev(Node* node, Node* before, Node* after);
    //void node_replace_next(Node* node, Node* before, Node* after);

    void to_dot(ostream& out);
    bool is_valid(void);

    void destroy_alignable_graph(void);
    gssw_graph* create_alignable_graph(
        int32_t match = 2,
        int32_t mismatch = 2,
        int32_t gap_open = 3,
        int32_t gap_extension = 1);
    void align(string& sequence);
    void align(Alignment& alignment);
    
    //void topological_sort(void); // badly needed

private:

    map<int64_t, gssw_node*> gssw_nodes;
    gssw_graph* _gssw_graph;
    int8_t* _gssw_nt_table;
    int8_t* _gssw_score_matrix;
    int32_t _gssw_match;
    int32_t _gssw_mismatch;
    int32_t _gssw_gap_open;
    int32_t _gssw_gap_extension;

};


} // end namespace vg

#endif