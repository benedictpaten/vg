#include "vg.h"

using namespace std;
using namespace google::protobuf;
using namespace vg;


//VariantGraph::VariantGraph(void) { };
// construct from protobufs
VariantGraph::VariantGraph(istream& in) {
    ParseFromIstream(&in);
    // populate by-id node index
    for (int64_t i = 0; i < nodes_size(); ++i) {
        Node* n = mutable_nodes(i);
        node_index[n] = i;
        node_by_id[n->id()] = n;
    }
    for (int64_t i = 0; i < edges_size(); ++i) {
        Edge* e = mutable_edges(i);
        edge_index[e] = i;
        edge_from_to[e->from()][e->to()] = e;
        edge_to_from[e->to()][e->from()] = e;
    }
}

VariantGraph::VariantGraph(vector<Node>& nodesv) {
    for (vector<Node>::iterator n = nodesv.begin(); n != nodesv.end(); ++n) {
        Node& node = *n;
        int64_t id = node.id();
        if (current_id() < id) {
            set_current_id(id);
        }
        Node* new_node = add_nodes(); // add it to the graph
        new_node->set_sequence(node.sequence());
        new_node->set_id(node.id());
        //Node& new_node = nodes(nodes_size()-1); // get a reference to it
        node_by_id[new_node->id()] = new_node; // and insert into our id lookup table
        node_index[new_node] = nodes_size()-1;
    }
}

// construct from VCF records
// --------------------------
// algorithm
// maintain a core reference path upon which we add new variants as they come
// addition procedure is the following
// find reference node overlapping our start position
// if it is already the end of a node, add the new node
// if it is not the end of a node, break it, insert edges from old->new
// go to end position of alt allele (could be the same position)
// if it already has a break, just point to the next node in line
// if it is not broken, break it and point to the next node
// add new node for alt alleles, connect to start and end node in reference path
// store the ref mapping as a property of the edges and nodes (this allows deletion edges and insertion subpaths)
//
VariantGraph::VariantGraph(vcf::VariantCallFile& variantCallFile, FastaReference& reference) {

    for (vector<string>::iterator r = reference.index->sequenceNames.begin();
         r != reference.index->sequenceNames.end(); ++r) {

        string& seqName = *r;
        map<long, Node*> reference_path;
        //map<long, set<Node*> > nodes; // for maintaining a reference-sorted graph
        string seq = reference.getSequence(seqName);

        Node* ref_node = create_node(seq);
        reference_path[0] = ref_node;

        // track the last nodes so that we can connect everything completely when variants occur in succession
        map<long, set<Node*> > nodes_by_end_position;

        variantCallFile.setRegion(seqName);
        vcf::Variant var(variantCallFile);

        while (variantCallFile.getNextVariant(var)) {

            int current_pos = (long int) var.position - 1;
            // decompose the alt
            bool flat_input_vcf = false; // hack
            map<string, vector<vcf::VariantAllele> > alternates = (flat_input_vcf ? var.flatAlternates() : var.parsedAlternates());
            for (map<string, vector<vcf::VariantAllele> >::iterator va = alternates.begin(); va !=alternates.end(); ++va) {
                vector<vcf::VariantAllele>& alleles = va->second;

                for (vector<vcf::VariantAllele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                    vcf::VariantAllele& allele = *a;

                    // reference alleles are provided naturally by the reference itself
                    if (allele.ref == allele.alt) {
                        continue;
                    }

                    long allele_start_pos = allele.position - 1;  // 0/1 based conversion... thanks vcflib!
                    long allele_end_pos = allele_start_pos + allele.ref.size();

                    if (allele_start_pos == 0) {
                        Node* root = create_node(""); // ensures that we can handle variation at first position (important when aligning)
                        reference_path[-1] = root;
                    }

                    Node* left_ref_node = NULL;
                    Node* middle_ref_node = NULL;
                    Node* right_ref_node = NULL;

                    // divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right) {
                    divide_path(reference_path,
                                allele_start_pos,
                                left_ref_node,
                                right_ref_node);

                    //cerr << "nodes: left: " << left_ref_node->id() << " right: " << right_ref_node->id() << endl;

                    // if the ref portion of the allele is not empty, then we need to make another cut
                    if (!allele.ref.empty()) {
                        divide_path(reference_path,
                                    allele_end_pos,
                                    middle_ref_node,
                                    right_ref_node);
                    }

                    Node* alt_node;
                    // create a new alt node and connect the pieces from before
                    if (!allele.alt.empty() && !allele.ref.empty()) {

                        alt_node = create_node(allele.alt);
                        //ref_map.add_node(alt_node, allele_start_pos, );
                        create_edge(left_ref_node, alt_node);
                        create_edge(alt_node, right_ref_node);

                        // XXXXXXXX middle is borked
                        // why do we have to force this edge back in?
                        // ... because it's not in the ref map?? (??)
                        create_edge(left_ref_node, middle_ref_node);

                        nodes_by_end_position[allele_end_pos].insert(alt_node);
                        nodes_by_end_position[allele_end_pos].insert(middle_ref_node);

                    } else if (!allele.alt.empty()) { // insertion

                        alt_node = create_node(allele.alt);
                        create_edge(left_ref_node, alt_node);
                        create_edge(alt_node, right_ref_node);
                        nodes_by_end_position[allele_end_pos].insert(alt_node);
                        nodes_by_end_position[allele_end_pos].insert(left_ref_node);

                    } else {// otherwise, we have a deletion

                        create_edge(left_ref_node, right_ref_node);
                        nodes_by_end_position[allele_end_pos].insert(left_ref_node);

                    }

                    if (allele_end_pos == seq.size()) {
                        // ensures that we can handle variation at first position (important when aligning)
                        Node* end = create_node("");
                        reference_path[allele_end_pos] = end;
                        create_edge(alt_node, end);
                        create_edge(middle_ref_node, end);
                    }

                    // if there are previous nodes, connect them
                    map<long, set<Node*> >::iterator ep = nodes_by_end_position.find(allele_start_pos);
                    if (ep != nodes_by_end_position.end()) {
                        set<Node*>& previous_nodes = ep->second;
                        for (set<Node*>::iterator n = previous_nodes.begin(); n != previous_nodes.end(); ++n) {
                            if (node_index.find(*n) != node_index.end()) {
                                if (middle_ref_node) {
                                    create_edge(*n, middle_ref_node);
                                }
                                create_edge(*n, alt_node);
                            }
                        }
                    }
                    // clean up previous
                    while (nodes_by_end_position.begin()->first < allele_start_pos) {
                        nodes_by_end_position.erase(nodes_by_end_position.begin()->first);
                    }

                    /*
                    if (!is_valid()) {
                        cerr << "graph is invalid after variant" << endl
                             << var << endl;
                        exit(1);
                    }
                    */
                }
            }
        }
    }
}

Edge* VariantGraph::create_edge(Node* from, Node* to) {
    return create_edge(from->id(), to->id());
}

Edge* VariantGraph::create_edge(int64_t from, int64_t to) {
    // prevent self-linking (violates DAG/partial ordering property)
    if (to == from) return NULL;
    // ensure the edge does not already exist
    Edge* edge = edge_from_to[from][to];
    if (edge) return edge;
    // if not, create it
    edge = add_edges();
    edge->set_from(from);
    edge->set_to(to);
    edge_from_to[from][to] = edge;
    edge_to_from[to][from] = edge;
    edge_index[edge] = edges_size()-1;
    return edge;
}

void VariantGraph::destroy_edge(Edge* edge) {
    //if (!is_valid()) cerr << "graph ain't valid" << endl;
    // erase from indexes
    edge_from_to[edge->from()].erase(edge->to());
    edge_to_from[edge->to()].erase(edge->from());

    // erase from edges by moving to end and dropping
    int lei = edges_size()-1;
    int tei = edge_index[edge];
    Edge* last = mutable_edges(lei);
    edge_index.erase(last);
    edge_index.erase(edge);

    // swap
    mutable_edges()->SwapElements(tei, lei);

    // point to new position
    Edge* nlast = mutable_edges(tei);

    // insert the new edge index position
    edge_index[nlast] = tei;

    // and fix edge indexes for moved edge object
    edge_from_to[nlast->from()][nlast->to()] = nlast;
    edge_to_from[nlast->to()][nlast->from()] = nlast;

    // drop the last position, erasing the node
    mutable_edges()->RemoveLast();

}

// use the VariantGraph class to generate ids
Node* VariantGraph::create_node(string seq) {
    // create the node
    Node* node = add_nodes();
    node->set_sequence(seq);
    node->set_id(current_id());
    set_current_id(current_id()+1);
    // copy it into the graph
    // and drop into our id index
    node_by_id[node->id()] = node;
    node_index[node] = nodes_size()-1;
    return node;
}

void VariantGraph::destroy_node(Node* node) {
    //if (!is_valid()) cerr << "graph is invalid before destroy_node" << endl;
    // remove edges associated with node
    set<Edge*> edges_to_destroy;
    map<int64_t, map<int64_t, Edge*> >::iterator e = edge_from_to.find(node->id());
    if (e != edge_from_to.end()) {
        for (map<int64_t, Edge*>::iterator f = e->second.begin();
             f != e->second.end(); ++f) {
            edges_to_destroy.insert(f->second);
        }
    }
    e = edge_to_from.find(node->id());
    if (e != edge_to_from.end()) {
        for (map<int64_t, Edge*>::iterator f = e->second.begin();
             f != e->second.end(); ++f) {
            edges_to_destroy.insert(f->second);
        }
    }
    for (set<Edge*>::iterator e = edges_to_destroy.begin();
         e != edges_to_destroy.end(); ++e) {
        destroy_edge(*e);
    }
    // assert cleanup
    edge_to_from.erase(node->id());
    edge_from_to.erase(node->id());

    // swap node with the last in nodes
    // call RemoveLast() to drop the node
    int lni = nodes_size()-1;
    int tni = node_index[node];
    Node* last = mutable_nodes(lni);
    mutable_nodes()->SwapElements(tni, lni);
    Node* nlast = mutable_nodes(tni);
    node_by_id[last->id()] = nlast;
    node_index.erase(last);
    node_index[nlast] = tni;
    node_by_id.erase(node->id());
    node_index.erase(node);
    mutable_nodes()->RemoveLast();
    //if (!is_valid()) cerr << "graph is invalid after destroy_node" << endl;
}

// utilities
void VariantGraph::divide_node(Node* node, int pos, Node*& left, Node*& right) {

    //cerr << "divide node " << node->id() << " @" << pos << endl;

    map<int64_t, map<int64_t, Edge*> >::iterator e;

    // make our left node
    left = create_node(node->sequence().substr(0,pos));

    // replace node connections to prev (left)
    e = edge_to_from.find(node->id());
    if (e != edge_to_from.end()) {
        for (map<int64_t, Edge*>::iterator p = e->second.begin();
             p != e->second.end(); ++p) {
            create_edge(p->first, left->id());
        }
    }

    // make our right node
    right = create_node(node->sequence().substr(pos,node->sequence().size()-1));

    // replace node connections to next (right)
    e = edge_from_to.find(node->id());
    if (e != edge_from_to.end()) {
        for (map<int64_t, Edge*>::iterator n = e->second.begin();
             n != e->second.end(); ++n) {
            create_edge(right->id(), n->first);
        }
    }

    // connect left to right
    create_edge(left, right);

    destroy_node(node);

}

// for dividing a path of nodes with an underlying coordinate system
void VariantGraph::divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right) {

    map<long, Node*>::iterator target = path.upper_bound(pos);
    --target; // we should now be pointing to the target ref node

    long node_pos = target->first;
    Node* old = target->second;
    
    // nothing to do
    if (node_pos == pos) {

        map<long, Node*>::iterator n = target; --n;
        left = n->second;
        right = target->second;

    } else {

        // divide the target node at our pos
        int diff = pos - node_pos;

        divide_node(old, diff, left, right);

        // left
        path[node_pos] = left;

        // right
        path[pos] = right;
    }
}

bool VariantGraph::is_valid(void) {
    for (int i = 0; i < nodes_size(); ++i) {
        Node* n = mutable_nodes(i);
    }
    for (int i = 0; i < edges_size(); ++i) {
        Edge* e = mutable_edges(i);
        int64_t f = e->from();
        int64_t t = e->to();
        if (node_by_id.find(f) == node_by_id.end()) {
            cerr << "graph invalid: edge index=" << i << " cannot find node (from) " << f << endl;
            return false;
        }
        if (node_by_id.find(t) == node_by_id.end()) {
            cerr << "graph invalid: edge index=" << i << " cannot find node (to) " << t << endl;
            return false;
        }
        if (edge_from_to.find(f) == edge_from_to.end()) {
            cerr << "graph invalid: edge index=" << i << " could not find entry in edges_from_to for node " << f << endl;
            return false;
        }
        if (edge_to_from.find(t) == edge_to_from.end()) {
            cerr << "graph invalid: edge index=" << i << " could not find entry in edges_to_from for node " << t << endl;
            return false;
        }
    }
    return true;
}

void VariantGraph::to_dot(ostream& out) {
    out << "digraph graphname {" << endl;
    out << "    node [shape=plaintext];" << endl;
    for (int i = 0; i < nodes_size(); ++i) {
        Node* n = mutable_nodes(i);
        out << "    " << n->id() << " [label=\"" << n->id() << ":" << n->sequence() << "\"];" << endl;
    }
    for (int i = 0; i < edges_size(); ++i) {
        Edge* e = mutable_edges(i);
        Node* p = node_by_id[e->from()];
        Node* n = node_by_id[e->to()];
        out << "    " << p->id() << " -> " << n->id() << ";" << endl;
    }
    out << "}" << endl;
}

void VariantGraph::destroy_alignable_graph(void) {
    if (_gssw_graph) gssw_graph_destroy(_gssw_graph);
    gssw_nodes.clear(); // these are freed via gssw_graph_destroy
    if (_gssw_nt_table) free(_gssw_nt_table);
    if (_gssw_score_matrix) free(_gssw_score_matrix);
}

gssw_graph* VariantGraph::create_alignable_graph(
    int32_t match,
    int32_t mismatch,
    int32_t gap_open,
    int32_t gap_extension
) {

    _gssw_match = match;
    _gssw_mismatch = mismatch;
    _gssw_gap_open = gap_open;
    _gssw_gap_extension = gap_extension;

    // these are used when setting up the nodes
    // they can be cleaned up via destroy_alignable_graph()
    _gssw_nt_table = gssw_create_nt_table();
	_gssw_score_matrix = gssw_create_score_matrix(_gssw_match, _gssw_mismatch);

    _gssw_graph = gssw_graph_create(nodes_size());

    for (int i = 0; i < nodes_size(); ++i) {
        Node* n = mutable_nodes(i);
        gssw_nodes[n->id()] = (gssw_node*)gssw_node_create(NULL, n->id(),
                                                           n->sequence().c_str(),
                                                           _gssw_nt_table,
                                                           _gssw_score_matrix);
    }

    for (int i = 0; i < edges_size(); ++i) {
        Edge* e = mutable_edges(i);
        gssw_nodes_add_edge(gssw_nodes[e->from()], gssw_nodes[e->to()]);
    }

    for (map<int64_t, gssw_node*>::iterator n = gssw_nodes.begin(); n != gssw_nodes.end(); ++n) {
        gssw_graph_add_node(_gssw_graph, n->second);
    }

}

void VariantGraph::align(string& sequence) {
    
    gssw_graph_fill(_gssw_graph, sequence.c_str(),
                    _gssw_nt_table, _gssw_score_matrix,
                    _gssw_gap_open, _gssw_gap_extension, 15, 2);

    gssw_graph_print_score_matrices(_gssw_graph, sequence.c_str(), sequence.size());
    gssw_graph_mapping* gm = gssw_graph_trace_back (_gssw_graph,
                                                    sequence.c_str(),
                                                    sequence.size(),
                                                    _gssw_match,
                                                    _gssw_mismatch,
                                                    _gssw_gap_open,
                                                    _gssw_gap_extension);

    gssw_print_graph_mapping(gm);
    gssw_graph_mapping_destroy(gm);

}

/*

	int32_t match = 2, mismatch = 2, gap_open = 3, gap_extension = 1;
    // from Mengyao's example about the importance of using all three matrices in traceback.
    // int32_t l, m, k, match = 2, mismatch = 1, gap_open = 2, gap_extension = 1;

    char *ref_seq_1 = argv[1];
    char *ref_seq_2 = argv[2];
    char *ref_seq_3 = argv[3];
    char *ref_seq_4 = argv[4];
    char *read_seq = argv[5];

    int8_t* nt_table = gssw_create_nt_table();
    
	// initialize scoring matrix for genome sequences
	//  A  C  G  T	N (or other ambiguous code)
	//  2 -2 -2 -2 	0	A
	// -2  2 -2 -2 	0	C
	// -2 -2  2 -2 	0	G
	// -2 -2 -2  2 	0	T
	//	0  0  0  0  0	N (or other ambiguous code)
	int8_t* mat = gssw_create_score_matrix(match, mismatch);

    gssw_node* nodes[4];
    nodes[0] = (gssw_node*)gssw_node_create("A", 1, ref_seq_1, nt_table, mat);
    nodes[1] = (gssw_node*)gssw_node_create("B", 2, ref_seq_2, nt_table, mat);
    nodes[2] = (gssw_node*)gssw_node_create("C", 3, ref_seq_3, nt_table, mat);
    nodes[3] = (gssw_node*)gssw_node_create("D", 4, ref_seq_4, nt_table, mat);

    // makes a diamond
    gssw_nodes_add_edge(nodes[0], nodes[1]);
    gssw_nodes_add_edge(nodes[0], nodes[2]);
    gssw_nodes_add_edge(nodes[1], nodes[3]);
    gssw_nodes_add_edge(nodes[2], nodes[3]);


    gssw_graph* graph = gssw_graph_create(4);
    //memcpy((void*)graph->nodes, (void*)nodes, 4*sizeof(gssw_node*));
    //graph->size = 4;
    gssw_graph_add_node(graph, nodes[0]);
    gssw_graph_add_node(graph, nodes[1]);
    gssw_graph_add_node(graph, nodes[2]);
    gssw_graph_add_node(graph, nodes[3]);

    gssw_graph_fill(graph, read_seq, nt_table, mat, gap_open, gap_extension, 15, 2);
    gssw_graph_print_score_matrices(graph, read_seq, strlen(read_seq));
    gssw_graph_mapping* gm = gssw_graph_trace_back (graph,
                                                    read_seq,
                                                    strlen(read_seq),
                                                    match,
                                                    mismatch,
                                                    gap_open,
                                                    gap_extension);

    gssw_print_graph_mapping(gm);
    gssw_graph_mapping_destroy(gm);
    // note that nodes which are referred to in this graph are destroyed as well
    gssw_graph_destroy(graph);

*/