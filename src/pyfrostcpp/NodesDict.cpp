#include "NodesDict.h"

namespace pyfrost {

void define_NodesDict(py::module& m) {
    auto py_NodesDict = py::class_<NodesDict>(m, "NodesDict")
        .def(py::init<PyfrostCCDBG&>())

        // Access nodes with [] operator overloading
        .def("__getitem__", py::overload_cast<char const*>(&NodesDict::findNode), py::is_operator())
        .def("__getitem__", py::overload_cast<Kmer const&>(&NodesDict::findNode), py::is_operator())

        .def("__contains__", py::overload_cast<Kmer const&>(&NodesDict::contains), py::is_operator())
        .def("__contains__", py::overload_cast<char const*>(&NodesDict::contains), py::is_operator())
        // Anything else than a k-mer is not present
        .def("__contains__", [] (NodesDict const& self, py::object const& o) { return false; })

        .def("__len__", &NodesDict::numNodes, py::is_operator())

        // Iterate over all nodes
        .def("__iter__", [](NodesDict const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    // NodesDict should act like a dict, so inherit from Mapping mixin
    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    auto Set = py::module::import("collections.abc").attr("Set");
    py_NodesDict.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(py_NodesDict.attr("__bases__"));
}

}
