#include <utility>
#include <vector>
#include <string>

#include "Pyfrost.h"
#include "NodeView.h"
#include "EdgeView.h"
#include "NodeBunchIter.h"
#include "AdjacencyProxy.h"

#ifndef PYFROST_BIFROSTGRAPH_H
#define PYFROST_BIFROSTGRAPH_H

namespace py = pybind11;

using std::vector;
using std::string;

namespace pyfrost {

class BifrostDiGraph {
public:
    friend class EdgeView;

    BifrostDiGraph() : nodes(dbg), edges(*this),
        succ(dbg, AdjacencyType::SUCCESSORS),
        pred(dbg, AdjacencyType::PREDECESSORS) {

        populateAttrs();
    }
    BifrostDiGraph(BifrostDiGraph const& o) : dbg(o.dbg), nodes(dbg), edges(*this),
        succ(dbg, AdjacencyType::SUCCESSORS),
        pred(dbg, AdjacencyType::PREDECESSORS),
        attr(o.attr) {

        populateAttrs();
    }

    BifrostDiGraph(BifrostDiGraph&& o) noexcept : dbg(std::move(o.dbg)), nodes(dbg), edges(*this),
        succ(dbg, AdjacencyType::SUCCESSORS),
        pred(dbg, AdjacencyType::PREDECESSORS),
        attr(std::move(o.attr)) {

        populateAttrs();
    }

    explicit BifrostDiGraph(PyfrostCCDBG&& o) noexcept : dbg(std::move(o)), nodes(dbg), edges(*this),
        succ(dbg, AdjacencyType::SUCCESSORS),
        pred(dbg, AdjacencyType::PREDECESSORS) {

        populateAttrs();
    }

    inline NodeView const& getNodeView() const {
        return nodes;
    }

    inline AdjacencyProxy const& getSuccessorsProxy() const {
        return succ;
    }

    inline AdjacencyProxy const& getPredecessorsProxy() const {
        return pred;
    }

    inline AdjacencyViewBase* getSuccessors(PyfrostColoredUMap const& unitig) {
        return succ.getView(unitig);
    }

    inline AdjacencyViewBase* getSuccessors(Kmer const& kmer) {
        return succ.getView(kmer);
    }

    inline AdjacencyViewBase* getSuccessors(char const* kmer) {
        return succ.getView(kmer);
    }

    inline AdjacencyViewBase* getPredecessors(PyfrostColoredUMap const& unitig) {
        return pred.getView(unitig);
    }

    inline AdjacencyViewBase* getPredecessors(Kmer const& kmer) {
        return pred.getView(kmer);
    }

    inline AdjacencyViewBase* getPredecessors(char const* kmer) {
        return pred.getView(kmer);
    }

    inline PyfrostColoredUMap findUnitig(Kmer const& kmer, bool extremities_only=false) {
        return dbg.find(kmer, extremities_only);
    }

    inline PyfrostColoredUMap findUnitig(char const* kmer, bool extremities_only=false) {
        return dbg.find(Kmer(kmer), extremities_only);
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

private:
    void populateAttrs() {
        attr["k"] = dbg.getK();
        attr["color_names"] = dbg.getData()->getColorNames();
    }

    PyfrostCCDBG dbg;
    NodeView nodes;
    EdgeView edges;
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
