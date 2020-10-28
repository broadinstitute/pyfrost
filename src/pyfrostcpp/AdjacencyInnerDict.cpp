#include "AdjacencyInnerDict.h"

namespace pyfrost {

/**
 * Returns the hardcoded metadata dict for an edge.
 */
py::dict makeEdgeDataDict(PyfrostCCDBG& dbg, Kmer const& kmer1, Kmer const& kmer2) {
    auto n1 = dbg.find(kmer1, true);
    auto n2 = dbg.find(kmer2, true);

    py::dict meta;
    if(n1.isEmpty || n2.isEmpty) {
        return meta;
    }

    meta["label"] = n2.getMappedHead().getChar(n2.getGraph()->getK() - 1);
    meta["orientation"] = py::make_tuple(
        n1.strand ? Strand::FORWARD : Strand::REVERSE,
        n2.strand ? Strand::FORWARD : Strand::REVERSE);

    return meta;
}

AdjacencyInnerDict::AdjacencyInnerDict(PyfrostCCDBG& _dbg, const Kmer& kmer, AdjacencyType _type)
    : dbg(_dbg), type(_type)
{
    node = getUnitig(kmer);
}


PyfrostColoredUMap AdjacencyInnerDict::getUnitig(Kmer const& kmer) {
    auto unitig = dbg.find(kmer, true).mappingToFullUnitig();

    if(unitig.isEmpty) {
        throw std::out_of_range("Node does not exist in the graph.");
    }

    return unitig;
}

AdjacencyInnerDict::iterator_type AdjacencyInnerDict::begin() const {
    return (type == AdjacencyType::SUCCESSORS
            ? iterator_type(&dbg, node.getSuccessors().begin(), false)
            : iterator_type(&dbg, node.getPredecessors().begin(), false));
}

AdjacencyInnerDict::iterator_type AdjacencyInnerDict::end() const {
    return (type == AdjacencyType::SUCCESSORS
            ? iterator_type(&dbg, node.getSuccessors().end(), false)
            : iterator_type(&dbg, node.getPredecessors().end(), false));
}

bool AdjacencyInnerDict::contains(Kmer const& o) const {
    // We have at most 4 neighbors so this should be quick enough, and could be considered constant time.
    return std::any_of(this->begin(), this->end(), [&] (Kmer const& neighbor) {
        return o == neighbor && !is_kmer_empty(neighbor);
    });
}

size_t AdjacencyInnerDict::numNeigbors() const {
    size_t num = 0;
    for(auto const& n : *this) {
        ++num;
    }

    return num;
}


py::dict AdjacencyInnerDict::getEdgeDict(Kmer const& neighbor) const {
    return (type == AdjacencyType::SUCCESSORS
            ? makeEdgeDataDict(dbg, node.getMappedHead(), neighbor)
            : makeEdgeDataDict(dbg, neighbor, node.getMappedHead()));
}

void define_AdjacencyInnerDict(py::module& m) {
    py::enum_<AdjacencyType>(m, "AdjacencyType")
        .value("SUCCESSORS", AdjacencyType::SUCCESSORS)
        .value("PREDECESSORS", AdjacencyType::PREDECESSORS);

    auto py_AdjacencyInnerDict = py::class_<AdjacencyInnerDict>(m, "AdjacencyInnerDict")
        .def("__iter__", [] (AdjacencyInnerDict const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        .def("__len__", &AdjacencyInnerDict::numNeigbors, py::is_operator())

        .def("__contains__", py::overload_cast<char const*>(&AdjacencyInnerDict::contains, py::const_),
            py::is_operator())
        .def("__contains__", py::overload_cast<Kmer const&>(&AdjacencyInnerDict::contains, py::const_),
            py::is_operator())
        .def("__contains__", [] (AdjacencyInnerDict const& self, py::object const& o) { return false; })

        .def("__getitem__", py::overload_cast<char const*>(&AdjacencyInnerDict::getEdgeDict, py::const_))
        .def("__getitem__", py::overload_cast<Kmer const&>(&AdjacencyInnerDict::getEdgeDict, py::const_));

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    auto Set = py::module::import("collections.abc").attr("Set");
    py_AdjacencyInnerDict.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(
        py_AdjacencyInnerDict.attr("__bases__"));

}


}
