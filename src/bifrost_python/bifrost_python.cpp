#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "Kmer.h"
#include "UnitigDataDict.h"
#include "UnitigMapping.h"
#include "NodeView.h"
#include "AdjacencyProxy.h"
#include "BifrostDiGraph.h"

namespace py = pybind11;


PYBIND11_MODULE(bifrost_python, m) {
    py::bind_vector<std::vector<std::string>>(m, "StringVector");

    m.doc() = R"doc(
        Python bindings for Bifrost
        ===========================

        This module provides a Python interface to the Bifrost colored compacted de Bruijn graph library, with a
        NetworkX compatible API.
    )doc";

    m.attr("default_k") = DEFAULT_K;

    pyfrost::define_Kmer(m);
    pyfrost::define_UnitigMapping(m);
    pyfrost::define_NodeView(m);
    pyfrost::define_AdjacencyProxy(m);
    pyfrost::define_BifrostDiGraph(m);

    m.def("reverse_complement", py::overload_cast<char const*>(&reverse_complement),
        "Return the reverse complement of a DNA string");
}
