#include <string>
#include <vector>

#include "Pyfrost.h"
#include "Kmer.h"
#include "KmerCounter.h"
#include "UnitigDataDict.h"
#include "UnitigDataProxy.h"
#include "UnitigMapping.h"
#include "UnitigColors.h"
#include "NodeView.h"
#include "AdjacencyProxy.h"
#include "DegreeView.h"
#include "BifrostDiGraph.h"

namespace py = pybind11;

std::ostream& operator<<(std::ostream& o, pyfrost::PyfrostColoredUMap const& u) {
    o << "<UnitigMap " << u.getMappedHead().toString() << ">";
    return o;
}

PYBIND11_MODULE(bifrost_python, m) {
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
    pyfrost::define_KmerCounter(m);
    pyfrost::define_UnitigColors(m);
    pyfrost::define_UnitigDataProxy(m);
    pyfrost::define_UnitigMapping(m);
    pyfrost::define_NodeView(m);
    pyfrost::define_EdgeView(m);
    pyfrost::define_AdjacencyProxy(m);
    pyfrost::define_DegreeView<void>(m, "WholeGraph");
    pyfrost::define_DegreeView<py::iterable>(m, "PyCollection");
    pyfrost::define_BifrostDiGraph(m);

    m.def("reverse_complement", py::overload_cast<char const*>(&reverse_complement),
        "Return the reverse complement of a DNA string");
}
