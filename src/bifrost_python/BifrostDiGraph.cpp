#include "BifrostDiGraph.h"

namespace pyfrost {

void populate_options(CCDBG_Build_opt &opt, py::kwargs const &kwargs) {
    if (kwargs.contains("k")) {
        opt.k = py::cast<int>(kwargs["k"]);
    }

    if (kwargs.contains("g")) {
        opt.g = py::cast<int>(kwargs["g"]);
    }

    if (kwargs.contains("threads")) {
        opt.nb_threads = py::cast<size_t>(kwargs["threads"]);
    }

    if (kwargs.contains("verbose")) {
        opt.verbose = py::cast<bool>(kwargs["verbose"]);
    }

    if (kwargs.contains("delete_isolated")) {
        opt.deleteIsolated = py::cast<bool>(kwargs["delete_isolated"]);
    }

    if (kwargs.contains("clip_tips")) {
        opt.clipTips = py::cast<bool>(kwargs["clip_tips"]);
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
    bool result = ccdbg.buildGraph(opt);

    if(!result) {
        throw std::runtime_error("Error building the graph.");
    }

    ccdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
    result = ccdbg.buildColors(opt);

    if(!result) {
        throw std::runtime_error("Error building coloring the graph.");
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

        .def("find", py::overload_cast<Kmer const&, bool>(&BifrostDiGraph::findUnitig), py::arg("kmer"),
             py::arg("extremities_only") = false)
        .def("find", py::overload_cast<char const*, bool>(&BifrostDiGraph::findUnitig), py::arg("kmer"),
             py::arg("extremities_only") = false)

            // neigbors, successors, and predececessors all should return an iterator
        .def("neighbors", [] (BifrostDiGraph& self, PyfrostColoredUMap const& node) {
            auto view = self.getSuccessors(node);
            return py::make_iterator<py::return_value_policy::copy>(view->begin(), view->end());
        })
        .def("neighbors", [] (BifrostDiGraph& self, Kmer const& node) {
            auto view = self.getSuccessors(node);
            return py::make_iterator<py::return_value_policy::copy>(view->begin(), view->end());
        })
        .def("neighbors", [] (BifrostDiGraph& self, char const* node) {
            auto view = self.getSuccessors(node);
            return py::make_iterator<py::return_value_policy::copy>(view->begin(), view->end());
        })

        .def("successors", [] (BifrostDiGraph& self, PyfrostColoredUMap const& node) {
            auto view = self.getSuccessors(node);
            return py::make_iterator<py::return_value_policy::copy>(view->begin(), view->end());
        })
        .def("successors", [] (BifrostDiGraph& self, Kmer const& node) {
            auto view = self.getSuccessors(node);
            return py::make_iterator<py::return_value_policy::copy>(view->begin(), view->end());
        })
        .def("successors", [] (BifrostDiGraph& self, char const* node) {
            auto view = self.getSuccessors(node);
            return py::make_iterator<py::return_value_policy::copy>(view->begin(), view->end());
        })

        .def("predecessors", [] (BifrostDiGraph& self, PyfrostColoredUMap const& node) {
            auto view = self.getPredecessors(node);
            return py::make_iterator<py::return_value_policy::copy>(view->begin(), view->end());
        })
        .def("predecessors", [] (BifrostDiGraph& self, Kmer const& node) {
            auto view = self.getPredecessors(node);
            return py::make_iterator<py::return_value_policy::copy>(view->begin(), view->end());
        })
        .def("predecessors", [] (BifrostDiGraph& self, char const* node) {
            auto view = self.getPredecessors(node);
            return py::make_iterator<py::return_value_policy::copy>(view->begin(), view->end());
        })

        .def("__len__", [] (BifrostDiGraph const& self) {
            return self.numNodes();
        })

        .def("__iter__", [] (BifrostDiGraph const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.getSuccessorsProxy().begin(),
                                                                    self.getSuccessorsProxy().end());
        })

        .def("nbunch_iter", [] (BifrostDiGraph const& self, const py::iterable& nbunch) {
            return py::make_iterator(NodeBunchIter(self, nbunch.begin()), NodeBunchIter(self, nbunch.end()));
        }, py::keep_alive<0, 1>())

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
                               "Get the graph attributes dictionary.")

        .def("is_directed", [] () { return true; })
        .def("is_multigraph", [] () { return false; })
        ;

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
