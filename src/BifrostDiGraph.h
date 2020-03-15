#include <vector>
#include <string>

#include <pybind11/pybind11.h>
#include <CompactedDBG.hpp>

#include "Pyfrost.h"
#include "NodeView.h"

#ifndef PYFROST_BIFROSTGRAPH_H
#define PYFROST_BIFROSTGRAPH_H

namespace py = pybind11;

using std::vector;
using std::string;

namespace pyfrost {

class BifrostDiGraph {
public:
    BifrostDiGraph() : nodes(dbg) { }
    BifrostDiGraph(BifrostDiGraph const& o) : dbg(o.dbg), nodes(dbg) { }
    BifrostDiGraph(BifrostDiGraph&& o) noexcept : dbg(std::move(o.dbg)), nodes(dbg) { }
    explicit BifrostDiGraph(PyfrostCCDBG&& dbg) noexcept : dbg(std::move(dbg)), nodes(dbg) { }

    inline NodeView& getNodeView() { return nodes; }

private:
    PyfrostCCDBG dbg;
    NodeView nodes;

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

        .def_property_readonly("nodes", &BifrostDiGraph::getNodeView, "Get nodes in the graph.");

    m.def("load", &load, "Load an existing colored Bifrost graph from a file.");
    m.def("build", &build, "Build a colored compacted Bifrost graph from references and sequencing data.");
    m.def("build_from_refs", &build_from_refs, "Build a colored compacted Bifrost graph from reference FASTA files.");
    m.def("build_from_samples", &build_from_samples, "Build a colored compacted Bifrost graph from sequencing sample files (FASTQ).");
}

}


#endif //PYFROST_BIFROSTGRAPH_H
