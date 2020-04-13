#include "EdgeView.h"
#include "BifrostDiGraph.h"

namespace pyfrost {


EdgeView::EdgeView(BifrostDiGraph& dbg) : dbg(dbg) { }

py::dict EdgeView::get(py::tuple const& edge) {
    py::dict metadata;

    PyfrostColoredUMap n1;
    PyfrostColoredUMap n2;
    if (py::isinstance<py::str>(edge[0])) {
        auto str = edge[0].cast<std::string>();
        n1 = dbg.dbg.find(Kmer(str.c_str()), true).mappingToFullUnitig();
    } else if (py::isinstance<Kmer>(edge[0])) {
        n1 = dbg.dbg.find(edge[0].cast<Kmer>(), true).mappingToFullUnitig();
    } else if (py::isinstance<PyfrostColoredUMap>(edge[0])) {
        n1 = edge[0].cast<PyfrostColoredUMap>();
    } else {
        throw py::type_error("Unsupported type given, only str, Kmer and UnitigMapping objects supported.");
    }

    if (py::isinstance<py::str>(edge[1])) {
        auto str = edge[1].cast<std::string>();
        n2 = dbg.dbg.find(Kmer(str.c_str()), true).mappingToFullUnitig();
    } else if (py::isinstance<Kmer>(edge[1])) {
        n2 = dbg.dbg.find(edge[1].cast<Kmer>(), true).mappingToFullUnitig();
    } else if (py::isinstance<PyfrostColoredUMap>(edge[1])) {
        n2 = edge[1].cast<PyfrostColoredUMap>();
    } else {
        throw py::type_error("Unsupported type given, only str, Kmer and UnitigMapping objects supported.");
    }

    if(n1.isEmpty || !n1.isFullMapping() || n2.isEmpty || !n2.isFullMapping()) {
        throw py::index_error("Edge does not exist in the graph");
    }

    return get(n1, n2);

}

py::dict EdgeView::get(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2) {
    py::dict meta;
    meta["label"] = n2.getMappedHead().getChar(dbg.dbg.getK() - 1);
    meta["orientation"] = py::make_tuple(
        n1.strand ? Strand::FORWARD : Strand::REVERSE,
        n2.strand ? Strand::FORWARD : Strand::REVERSE
    );

    return meta;
}

OutEdgeIterator<PyfrostCCDBG::iterator> EdgeView::begin() const {
    return OutEdgeIterator<PyfrostCCDBG::iterator>(dbg.getNodeView().begin());
}

OutEdgeIterator<PyfrostCCDBG::iterator> EdgeView::end() const {
    return OutEdgeIterator<PyfrostCCDBG::iterator>(dbg.getNodeView().end());
}

size_t EdgeView::size() const {
    size_t num_edges = 0;
    for(auto const& um : dbg.getNodeView()) {
        num_edges += um.getSuccessors().cardinality();
    }

    return num_edges;
}


void define_EdgeView(py::module& m) {
    auto py_EdgeView = py::class_<EdgeView>(m, "EdgeView")
        .def("__getitem__", py::overload_cast<py::tuple const&>(&EdgeView::get), py::is_operator())
        .def("__len__", &EdgeView::size)
        .def("__iter__", [] (EdgeView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    py_EdgeView.attr("__bases__") = py::make_tuple(Mapping).attr("__add__")(py_EdgeView.attr("__bases__"));

}

}
