#include "EdgeView.h"
#include "BifrostDiGraph.h"

namespace pyfrost {

/**
 * Returns the hardcoded metadata dict for an edge.
 */
py::dict makeEdgeDataDict(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2) {
    py::dict meta;
    if(n1.isEmpty || n2.isEmpty) {
        return meta;
    }

    meta["label"] = n2.getMappedHead().getChar(n2.getGraph()->getK() - 1);
    meta["orientation"] = py::make_tuple(
        n1.strand ? Strand::FORWARD : Strand::REVERSE,
        n2.strand ? Strand::FORWARD : Strand::REVERSE
    );

    return meta;
}

PyfrostColoredUMap findNode(PyfrostCCDBG& dbg, py::object const& n) {
    PyfrostColoredUMap node;

    if (py::isinstance<py::str>(n)) {
        auto str = n.cast<std::string>();
        node = dbg.find(Kmer(str.c_str()), true).mappingToFullUnitig();
    } else if (py::isinstance<Kmer>(n)) {
        node = dbg.find(n.cast<Kmer>(), true).mappingToFullUnitig();
    } else if (py::isinstance<PyfrostColoredUMap>(n)) {
        node = n.cast<PyfrostColoredUMap>();
    } else {
        throw py::type_error("Unsupported type given, only str, Kmer and UnitigMapping objects supported.");
    }

    if(node.isEmpty || !node.isFullMapping()) {
        throw py::index_error("Node does not exists in the graph");
    }

    return node;
}

bool isSuccessor(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2) {
    for(auto const& succ : n1.getSuccessors()) {
        if(succ == n2) {
            return true;
        }
    }

    return false;
}

bool isPredecessor(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2) {
    for(auto const& pred : n1.getPredecessors()) {
        if(pred == n2) {
            return true;
        }
    }

    return false;
}

//
// OutEdgeView
// -----------
//

OutEdgeView::OutEdgeView(BifrostDiGraph& dbg) : dbg(dbg) { }

py::dict OutEdgeView::get(py::tuple const& edge) {
    py::dict metadata;

    PyfrostColoredUMap n1 = findNode(dbg.dbg, edge[0]);
    PyfrostColoredUMap n2 = findNode(dbg.dbg, edge[1]);

    return get(n1, n2);
}

py::dict OutEdgeView::get(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2) {
    if(!isSuccessor(n1, n2)) {
        throw py::index_error("Target node is not a successor of source node.");
    }

    return makeEdgeDataDict(n1, n2);
}

EdgeIterator<PyfrostCCDBG::iterator> OutEdgeView::begin() const {
    return {dbg.getNodeView().begin(), dbg.getNodeView().end()};
}

EdgeIterator<PyfrostCCDBG::iterator> OutEdgeView::end() const {
    return {};
}

bool OutEdgeView::contains(const py::tuple &edge) {
    auto n1 = findNode(dbg.dbg, edge[0]);
    auto n2 = findNode(dbg.dbg, edge[1]);

    return isSuccessor(n1, n2);
}

size_t OutEdgeView::size() const {
    size_t num_edges = 0;
    for(auto const& um : dbg.getNodeView()) {
        num_edges += um.getSuccessors().cardinality();
    }

    return num_edges;
}

//
// InEdgeView
// -----------
//

InEdgeView::InEdgeView(BifrostDiGraph& dbg) : dbg(dbg) { }

py::dict InEdgeView::get(py::tuple const& edge) {
    PyfrostColoredUMap n1 = findNode(dbg.dbg, edge[0]);
    PyfrostColoredUMap n2 = findNode(dbg.dbg, edge[1]);

    return get(n1, n2);
}

py::dict InEdgeView::get(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2) {
    if(!isSuccessor(n1, n2)) {
        throw py::index_error("Edge not does exists");
    }

    return makeEdgeDataDict(n1, n2);
}

EdgeIterator<PyfrostCCDBG::iterator, false> InEdgeView::begin() const {
    return {dbg.getNodeView().begin(), dbg.getNodeView().end()};
}

EdgeIterator<PyfrostCCDBG::iterator, false> InEdgeView::end() const {
    return {};
}

bool InEdgeView::contains(const py::tuple &edge) {
    auto n1 = findNode(dbg.dbg, edge[0]);
    auto n2 = findNode(dbg.dbg, edge[1]);

    return isSuccessor(n1, n2);
}

size_t InEdgeView::size() const {
    size_t num_edges = 0;
    for(auto const& um : dbg.getNodeView()) {
        num_edges += um.getPredecessors().cardinality();
    }

    return num_edges;
}


void define_EdgeView(py::module& m) {
    auto py_OutEdgeView = py::class_<OutEdgeView>(m, "OutEdgeView")
        .def("__getitem__", py::overload_cast<py::tuple const&>(&OutEdgeView::get), py::is_operator())
        .def("__contains__", &OutEdgeView::contains, py::is_operator())
        .def("__len__", &OutEdgeView::size)
        .def("__iter__", [] (OutEdgeView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    auto Set = py::module::import("collections.abc").attr("Set");
    py_OutEdgeView.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(py_OutEdgeView.attr("__bases__"));

    auto py_InEdgeView = py::class_<InEdgeView>(m, "InEdgeView")
        .def("__getitem__", py::overload_cast<py::tuple const&>(&InEdgeView::get), py::is_operator())
        .def("__contains__", &InEdgeView::contains, py::is_operator())
        .def("__len__", &InEdgeView::size)
        .def("__iter__", [] (InEdgeView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    py_InEdgeView.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(py_InEdgeView.attr("__bases__"));

}

}
