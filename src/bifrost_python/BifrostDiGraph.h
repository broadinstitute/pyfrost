#include <utility>
#include <vector>
#include <string>

#include "Pyfrost.h"
#include "NodeView.h"
#include "EdgeView.h"
#include "AdjacencyProxy.h"

#ifndef PYFROST_BIFROSTGRAPH_H
#define PYFROST_BIFROSTGRAPH_H

namespace py = pybind11;

using std::vector;
using std::string;

namespace pyfrost {

class BifrostDiGraph {
public:
    BifrostDiGraph() : nodes(dbg), out_edges(dbg), in_edges(dbg),
                       succ(dbg, AdjacencyType::SUCCESSORS),
                       pred(dbg, AdjacencyType::PREDECESSORS) {

        populateAttrs();
    }
    BifrostDiGraph(BifrostDiGraph const& o) : dbg(o.dbg), nodes(dbg), out_edges(dbg), in_edges(dbg),
                                              succ(dbg, AdjacencyType::SUCCESSORS),
                                              pred(dbg, AdjacencyType::PREDECESSORS),
                                              attr(o.attr) {

        populateAttrs();
    }

    BifrostDiGraph(BifrostDiGraph&& o) noexcept : dbg(std::move(o.dbg)), nodes(dbg),
                                                  out_edges(dbg), in_edges(dbg),
                                                  succ(dbg, AdjacencyType::SUCCESSORS),
                                                  pred(dbg, AdjacencyType::PREDECESSORS),
                                                  attr(std::move(o.attr)) {

        populateAttrs();
    }

    explicit BifrostDiGraph(PyfrostCCDBG&& o) noexcept : dbg(std::move(o)), nodes(dbg),
                                                         out_edges(dbg), in_edges(dbg),
                                                         succ(dbg, AdjacencyType::SUCCESSORS),
                                                         pred(dbg, AdjacencyType::PREDECESSORS) {

        populateAttrs();
    }

    inline bool remove(Kmer const& kmer) {
        auto unitig = dbg.find(kmer, true);
        if(unitig.isEmpty) {
            throw py::value_error("Kmer is not a unitig.");
        }

        return dbg.remove(unitig);
    }

    inline bool remove(char const* kmer) {
        return remove(Kmer(kmer));
    }

    inline NodeView const& getNodeView() const {
        return nodes;
    }

    inline EdgeView<true> const& getOutEdgeView() const {
        return out_edges;
    }

    inline EdgeView<false> const& getInEdgeView() const {
        return in_edges;
    }

    inline AdjacencyProxy const& getSuccessorsProxy() const {
        return succ;
    }

    inline AdjacencyProxy const& getPredecessorsProxy() const {
        return pred;
    }

    inline SuccessorView getSuccessors(Kmer const& kmer) {
        return SuccessorView(dbg, kmer);
    }

    inline SuccessorView getSuccessors(char const* kmer) {
        return getSuccessors(Kmer(kmer));
    }

    inline PredecessorView getPredecessors(Kmer const& kmer) {
        return PredecessorView(dbg, kmer);
    }

    inline PredecessorView getPredecessors(char const* kmer) {
        return getPredecessors(Kmer(kmer));
    }

    inline UnitigDataProxy findUnitig(Kmer const& kmer, bool extremities_only=false) {
        auto unitig = dbg.find(kmer, extremities_only);

        if(unitig.isEmpty) {
            throw py::key_error("K-mer not found in the graph.");
        }

        return UnitigDataProxy(unitig);
    }

    inline UnitigDataProxy findUnitig(char const* kmer, bool extremities_only=false) {
        if(std::strlen(kmer) != Kmer::k) {
            throw py::value_error("Given k-mer is is not of length k!");
        }

        return findUnitig(Kmer(kmer), extremities_only);
    }

    inline bool contains(Kmer const& kmer) {
        return !dbg.find(kmer, true).isEmpty;
    }

    inline bool contains(char const* kmer) {
        return contains(Kmer(kmer));
    }

    inline size_t numNodes() const {
        return dbg.size();
    }

    inline py::dict& getGraphAttributes() {
        return attr;
    }

    bool operator==(BifrostDiGraph const& o) const {
        return dbg == o.dbg;
    }

    bool operator!=(BifrostDiGraph const& o) const {
        return !operator==(o);
    }

    inline PyfrostCCDBG& getGraph() {
        return dbg;
    }

private:
    void populateAttrs() {
        attr["k"] = dbg.getK();
        attr["g"] = dbg.getG();
        attr["color_names"] = dbg.getData()->getColorNames();
    }

    PyfrostCCDBG dbg;
    NodeView nodes;
    EdgeView<true> out_edges;
    EdgeView<false> in_edges;
    AdjacencyProxy succ;
    AdjacencyProxy pred;
    py::dict attr;

};


void populate_options(CCDBG_Build_opt& opt, py::kwargs const& kwargs);

BifrostDiGraph build(py::list const& input_seq_files, py::list const& input_ref_files,
                     py::kwargs const& kwargs);

BifrostDiGraph build_from_refs(py::list const& input_ref_files, py::kwargs const& kwargs);

BifrostDiGraph build_from_samples(py::list const& input_seq_files, py::kwargs const& kwargs);


BifrostDiGraph load(char const* input_graph_file, char const* input_color_file, py::kwargs const& kwargs);

void define_BifrostDiGraph(py::module& m);

}

#endif //PYFROST_BIFROSTGRAPH_H
