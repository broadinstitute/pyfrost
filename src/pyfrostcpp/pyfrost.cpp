#include <string>
#include <vector>

#include "pyfrost.h"
#include "Kmer.h"
#include "KmerCounter.h"
#include "UnitigDataDict.h"
#include "NodeDataDict.h"
#include "UnitigMapping.h"
#include "UnitigColors.h"
#include "NodesDict.h"
#include "AdjacencyOuterDict.h"
#include "JunctionTree.h"
#include "LinkDB.h"
#include "MemLinkDB.h"

namespace py = pybind11;

std::ostream& operator<<(std::ostream& o, pyfrost::PyfrostColoredUMap const& u) {
    o << "<UnitigMap " << u.getMappedHead().toString() << ">";
    return o;
}

namespace pyfrost {

NodeDataDict findUnitig(PyfrostCCDBG& g, Kmer const& kmer, bool extremities_only=false) {
    auto unitig = g.find(kmer, extremities_only);

    if(unitig.isEmpty) {
        throw py::value_error("Could not find k-mer in the graph.");
    }

    return NodeDataDict(unitig);
}

NodeDataDict findUnitig(PyfrostCCDBG& g, char const* kmer, bool extremities_only=false) {
    return findUnitig(g, Kmer(kmer), extremities_only);
}

void removeUnitig(PyfrostCCDBG& g, Kmer const& kmer) {
    auto unitig = g.find(kmer, true);

    if(unitig.isEmpty) {
        throw py::key_error("Node to remove doesn't exists");
    }

    g.remove(unitig);
}

void removeUnitig(PyfrostCCDBG& g, char const* kmer) {
    return removeUnitig(g, Kmer(kmer));
}

void define_PyfrostCCDBG(py::module& m) {
    py::class_<PyfrostCCDBG>(m, "PyfrostCCDBG")
        .def("get_k", [] (PyfrostCCDBG const& self) {
            return self.getK();
        })
        .def("get_g", [] (PyfrostCCDBG const& self) {
            return self.getG();
        })
        .def("color_names", [] (PyfrostCCDBG const& self) {
            return self.getData()->getColorNames();
        })
        .def("find", py::overload_cast<PyfrostCCDBG&, Kmer const&, bool>(&findUnitig),
            py::arg("kmer"), py::arg("extremities_only") = false)
        .def("find", py::overload_cast<PyfrostCCDBG&, char const*, bool>(&findUnitig),
            py::arg("kmer"), py::arg("extremities_only") = false)

        .def("remove", py::overload_cast<PyfrostCCDBG&, Kmer const&>(&removeUnitig))
        .def("remove", py::overload_cast<PyfrostCCDBG&, char const*>(&removeUnitig))
        ;
}

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


PyfrostCCDBG build(py::list const& input_ref_files, py::list const& input_seq_files,
                   py::kwargs const& kwargs) {
    CCDBG_Build_opt opt;

    if(!input_ref_files.empty()) {
        for(auto const& item : input_ref_files) {
            opt.filename_ref_in.emplace_back(py::cast<string>(item));
        }
    }

    if(!input_seq_files.empty()) {
        for(auto const& item : input_seq_files) {
            opt.filename_seq_in.emplace_back(py::cast<string>(item));
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

    return ccdbg;
}

PyfrostCCDBG load(char const* input_graph_file, char const* input_color_file, py::kwargs const& kwargs) {
    CCDBG_Build_opt opt;
    opt.filename_graph_in = input_graph_file;
    opt.filename_colors_in = input_color_file;

    populate_options(opt, kwargs);

    PyfrostCCDBG ccdbg(opt.k, opt.g);
    bool result = ccdbg.read(opt.filename_graph_in, opt.filename_colors_in, opt.nb_threads, opt.verbose);

    if(!result) {
        throw std::runtime_error("Error reading the graph.");
    }

    return ccdbg;
}

void dump(PyfrostCCDBG& g, string const& fname_prefix, size_t num_threads=2) {
    g.write(fname_prefix, num_threads);
    // TODO: write additional GFA tags
}

}

PYBIND11_MODULE(pyfrostcpp, m) {
    m.doc() = R"doc(
        Python bindings for Bifrost
        ===========================

        This module provides a Python interface to the Bifrost colored compacted de Bruijn graph library, with a
        NetworkX compatible API.
    )doc";

    m.attr("default_k") = DEFAULT_K;

    py::bind_vector<std::vector<std::string>>(m, "StringVector");
    py::implicitly_convertible<py::list, std::vector<std::string>>();

    pyfrost::define_Kmer(m);
    pyfrost::define_PyfrostCCDBG(m);
    pyfrost::define_KmerCounter(m);
    pyfrost::define_UnitigColors(m);
    pyfrost::define_NodeDataDict(m);
    pyfrost::define_UnitigMapping(m);
    pyfrost::define_NodesDict(m);
    pyfrost::define_AdjacencyInnerDict(m);
    pyfrost::define_AdjacencyOuterDict(m);
    pyfrost::define_JunctionTreeNode(m);
    pyfrost::define_LinkDB(m);
    pyfrost::define_MemLinkDB(m);

    m.def("load", &pyfrost::load,
          "Load an existing colored Bifrost graph from a file.");
    m.def("build", &pyfrost::build,
          "Build a colored compacted Bifrost graph from references and sequencing data.");
    m.def("dump", &pyfrost::dump, py::arg("g"), py::arg("fname_prefix"), py::arg("num_threads") = 2,
          "Save graph to file.");

    m.def("reverse_complement", py::overload_cast<char const*>(&reverse_complement),
        "Return the reverse complement of a DNA string");

    m.def("k_g", [] () {
        return py::make_tuple(Kmer::k, Minimizer::g);
    });
}
