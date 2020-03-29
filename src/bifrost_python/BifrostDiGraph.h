#include <vector>
#include <string>

#include <pybind11/pybind11.h>
#include <CompactedDBG.hpp>

#include "Pyfrost.h"
#include "NodeView.h"
#include "AdjacencyProxy.h"

#ifndef PYFROST_BIFROSTGRAPH_H
#define PYFROST_BIFROSTGRAPH_H

namespace py = pybind11;

using std::vector;
using std::string;

namespace pyfrost {

class BifrostDiGraph {
public:
    BifrostDiGraph() : nodes(dbg),
        succ(dbg, AdjacencyType::SUCCESSORS),
        pred(dbg, AdjacencyType::PREDECESSORS) { }
    BifrostDiGraph(BifrostDiGraph const& o) : dbg(o.dbg), nodes(dbg),
        succ(dbg, AdjacencyType::SUCCESSORS),
        pred(dbg, AdjacencyType::PREDECESSORS),
        attr(o.attr) { }
    BifrostDiGraph(BifrostDiGraph&& o) noexcept : dbg(std::move(o.dbg)), nodes(dbg),
        succ(dbg, AdjacencyType::SUCCESSORS),
        pred(dbg, AdjacencyType::PREDECESSORS),
        attr(std::move(o.attr)) { }

    explicit BifrostDiGraph(PyfrostCCDBG&& o) noexcept : dbg(std::move(o)), nodes(dbg),
        succ(dbg, AdjacencyType::SUCCESSORS),
        pred(dbg, AdjacencyType::PREDECESSORS) { }

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

    inline size_t numNodes() const {
        return dbg.size();
    }

    inline py::dict& getGraphAttributes() {
        return attr;
    }

private:
    PyfrostCCDBG dbg;
    NodeView nodes;
    AdjacencyProxy succ;
    AdjacencyProxy pred;
    py::dict attr;

};

void populate_options(CCDBG_Build_opt& opt, py::kwargs const& kwargs, bool build_input_files=false) {
    if(kwargs.contains("k")) {
        opt.k = py::cast<int>(kwargs["k"]);
    }

    if(kwargs.contains("g")) {
        opt.g = py::cast<int>(kwargs["g"]);
    }

    if(kwargs.contains("threads")) {
        opt.nb_threads = py::cast<size_t>(kwargs["threads"]);
    }

    if(kwargs.contains("verbose")) {
        opt.verbose = py::cast<bool>(kwargs["bool"]);
    }
}

BifrostDiGraph build(py::list const& input_seq_files, py::list const& input_ref_files,
                     py::kwargs const& kwargs) {
    CCDBG_Build_opt opt;
    if(!input_seq_files.empty()) {
        for(auto const& item : input_seq_files) {
            opt.filename_seq_in.emplace_back(py::cast<string>(item));
        }
    }

    if(!input_ref_files.empty()) {
        for(auto const& item : input_ref_files) {
            opt.filename_ref_in.emplace_back(py::cast<string>(item));
        }
    }

    populate_options(opt, kwargs);

    PyfrostCCDBG ccdbg(opt.k, opt.g);
    bool result = ccdbg.build(opt);

    if(!result) {
        throw std::runtime_error("Error building the graph.");
    }

    return BifrostDiGraph(std::move(ccdbg));
}

BifrostDiGraph build_from_refs(py::list const& input_ref_files, py::kwargs const& kwargs) {

    return build(py::list(), input_ref_files, kwargs);
}

BifrostDiGraph build_from_samples(py::list const& input_seq_files, py::kwargs const& kwargs) {

    return build(input_seq_files, py::list(), kwargs);
}


BifrostDiGraph load(char const* input_graph_file, char const* input_color_file, py::kwargs const& kwargs) {
    CCDBG_Build_opt opt;
    opt.filename_graph_in = input_graph_file;
    opt.filename_colors_in = input_color_file;

    populate_options(opt, kwargs);

    PyfrostCCDBG ccdbg(opt.k, opt.g);
    bool result = ccdbg.read(opt.filename_graph_in, opt.filename_colors_in, opt.nb_threads, opt.verbose);

    if(!result) {
        throw std::runtime_error("Error reading the graph.");
    }

    return BifrostDiGraph(std::move(ccdbg));
}

void define_BifrostDiGraph(py::module& m) {
    py::class_<BifrostDiGraph>(m, "BifrostDiGraph")
        .def(py::init<>())
        .def(py::init<BifrostDiGraph const&>())

        .def("__getitem__", py::overload_cast<PyfrostColoredUMap const&>(&BifrostDiGraph::getSuccessors),
            "Access a node")
        .def("__getitem__", py::overload_cast<Kmer const&>(&BifrostDiGraph::getSuccessors),
             "Access a node")
        .def("__getitem__", py::overload_cast<char const*>(&BifrostDiGraph::getSuccessors),
             "Access a node")

        .def("neighbors", py::overload_cast<PyfrostColoredUMap const&>(&BifrostDiGraph::getSuccessors),
             "Access a node")
        .def("neighbors", py::overload_cast<Kmer const&>(&BifrostDiGraph::getSuccessors),
             "Access a node")
        .def("neighbors", py::overload_cast<char const*>(&BifrostDiGraph::getSuccessors),
             "Access a node")

        .def("successors", py::overload_cast<PyfrostColoredUMap const&>(&BifrostDiGraph::getSuccessors),
             "Access a node")
        .def("successors", py::overload_cast<Kmer const&>(&BifrostDiGraph::getSuccessors),
             "Access a node")
        .def("successors", py::overload_cast<char const*>(&BifrostDiGraph::getSuccessors),
             "Access a node")

        .def("predecessors", py::overload_cast<PyfrostColoredUMap const&>(&BifrostDiGraph::getPredecessors),
             "Access a node")
        .def("predecessors", py::overload_cast<Kmer const&>(&BifrostDiGraph::getPredecessors),
             "Access a node")
        .def("predecessors", py::overload_cast<char const*>(&BifrostDiGraph::getPredecessors),
             "Access a node")

        .def("__len__", [] (BifrostDiGraph const& self) {
            return self.numNodes();
        })

        .def("__iter__", [] (BifrostDiGraph const& self) {
            return py::make_iterator(self.getSuccessorsProxy().begin(), self.getSuccessorsProxy().end());
        })

        .def("number_of_nodes", &BifrostDiGraph::numNodes, "Get the number of nodes.")

        .def_property_readonly("nodes", &BifrostDiGraph::getNodeView,
            "Get nodes in the graph.")
        .def_property_readonly("succ", &BifrostDiGraph::getSuccessorsProxy,
            "Get successors keyed by node.")
        .def_property_readonly("adj", &BifrostDiGraph::getSuccessorsProxy,
            "Get successors keyed by node.")
        .def_property_readonly("pred", &BifrostDiGraph::getPredecessorsProxy,
            "Get predecessors keyed by node.")
        .def_property_readonly("graph", &BifrostDiGraph::getGraphAttributes,
            "Get the graph attributes dictionary.");

    m.def("load", &load, py::return_value_policy::move,
        "Load an existing colored Bifrost graph from a file.");
    m.def("build", &build, py::return_value_policy::move,
        "Build a colored compacted Bifrost graph from references and sequencing data.");
    m.def("build_from_refs", &build_from_refs, py::return_value_policy::move,
        "Build a colored compacted Bifrost graph from reference FASTA files.");
    m.def("build_from_samples", &build_from_samples, py::return_value_policy::move,
        "Build a colored compacted Bifrost graph from sequencing sample files (FASTQ).");
}

}


#endif //PYFROST_BIFROSTGRAPH_H
