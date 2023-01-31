#include <string>
#include <vector>

#include "pyfrost.h"
#include "Kmer.h"
#include "Minimizers.h"
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
#include "LinkAnnotator.h"
#include "Neighbors.h"

namespace py = pybind11;

std::ostream& operator<<(std::ostream& o, pyfrost::PyfrostColoredUMap const& u) {
    o << "<UnitigMap " << u.getMappedHead().toString() << ">";
    return o;
}

namespace pyfrost {

void setKG(size_t k, size_t g)
{
    if(k <= 2) {
        throw std::out_of_range("k-mer size needs to be at least 3");
    }

    if(k >= MAX_KMER_SIZE) {
        const uint16_t max_kmer_size = MAX_KMER_SIZE - 1;
        std::stringstream error_msg;
        error_msg << "K-mer size is too big! Max k-mer size: " << max_kmer_size;
        throw std::out_of_range(error_msg.str());
    }

    if(g > (k - 2)) {
        throw std::out_of_range("Minimizer length cannot exceed k-2");
    }

    if(g == 0) {
        if(k >= 15) {
            g = k - DEFAULT_G_DEC1;
        } else if(k >= 7) {
            g = k - DEFAULT_G_DEC2;
        } else {
            g = k - 2;
        }
    }

    g = std::min(MAX_GMER_SIZE, (int) g);

    if(Kmer::k > 0 && Kmer::k != k) {
        std::cerr << "WARNING: setting new k-mer size! old: " << Kmer::k << " => new: " << k << std::endl;
    }
    Kmer::set_k(k);

    if(Minimizer::g > 0 && Minimizer::g != g) {
        std::cerr << "WARNING: setting new minimizer size! old: " << Minimizer::g << " => new: " << g << std::endl;
    }
    Minimizer::set_g(g);
}

py::object findUnitig(PyfrostCCDBG& g, Kmer const& kmer, bool extremities_only=false) {
    auto unitig = g.find(kmer, extremities_only);

    if(unitig.isEmpty) {
        return py::none();
    }

    return py::cast(NodeDataDict(unitig));
}

py::object findUnitig(PyfrostCCDBG& g, char const* kmer, bool extremities_only=false) {
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

        .def("color_restricted_successors", [] (PyfrostCCDBG& self, Kmer const& node,
                                                unordered_set<size_t> const& allowed_colors) {
            auto unitig = self.find(node, true);

            return colorRestrictedSuccessors(unitig, allowed_colors);
        })
        .def("color_restricted_predecessors", [] (PyfrostCCDBG& self, Kmer const& node,
                unordered_set<size_t> const& allowed_colors) {
            auto unitig = self.find(node, true);

            return colorRestrictedPredecessors(unitig, allowed_colors);
        })
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
        throw std::runtime_error("Error building coloring of the graph.");
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
    py::bind_map<std::unordered_map<Kmer, size_t>>(m, "KmerSizeTMap");
    py::implicitly_convertible<py::list, std::vector<std::string>>();

    pyfrost::define_Kmer(m);
    pyfrost::define_KmerCounter(m);

    pyfrost::define_Minimizer(m);
    pyfrost::define_MinHashIterator(m);
    pyfrost::define_MinHashResult(m);

    pyfrost::define_PyfrostCCDBG(m);
    pyfrost::define_UnitigColors(m);
    pyfrost::define_NodeDataDict(m);
    pyfrost::define_UnitigMapping(m);
    pyfrost::define_NodesDict(m);
    pyfrost::define_AdjacencyInnerDict(m);
    pyfrost::define_AdjacencyOuterDict(m);

    pyfrost::define_JunctionTreeNode(m);
    pyfrost::define_LinkDB(m);
    pyfrost::define_MemLinkDB<pyfrost::JunctionTreeNode>(m, "MemLinkDB");
    pyfrost::define_MemLinkDB<pyfrost::JunctionTreeNodeWithCov>(m, "MemLinkDBWithCov");
    pyfrost::define_LinkAnnotator(m);
    pyfrost::define_MappingResult(m);

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
