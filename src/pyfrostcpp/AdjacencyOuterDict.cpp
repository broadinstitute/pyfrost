#include "AdjacencyOuterDict.h"

namespace pyfrost {

void define_AdjacencyOuterDict(py::module& m) {

    auto py_AdjacencyOuterDict = py::class_<AdjacencyOuterDict>(m, "AdjacencyOuterDict")
        .def(py::init<PyfrostCCDBG&, AdjacencyType>())

        .def("__len__", &AdjacencyOuterDict::numNodes, py::is_operator())

        .def("__contains__", py::overload_cast<Kmer const&>(&AdjacencyOuterDict::contains))
        .def("__contains__", py::overload_cast<char const*>(&AdjacencyOuterDict::contains))
        .def("__contains__", [] (AdjacencyOuterDict const& self, py::object const& o) { return false; })

        .def("__getitem__", py::overload_cast<Kmer const&>(&AdjacencyOuterDict::getInnerDict), py::is_operator(),
             py::keep_alive<0, 1>())
        .def("__getitem__", py::overload_cast<char const*>(&AdjacencyOuterDict::getInnerDict), py::is_operator(),
             py::keep_alive<0, 1>())

        .def("__iter__", [](AdjacencyOuterDict const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.begin(), self.end());
        }, py::keep_alive<0, 1>())
        .def("iter_no_rev_compl", [](AdjacencyOuterDict const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.begin_no_rc(), self.end_no_rc());
        }, py::keep_alive<0, 1>());

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    auto Set = py::module::import("collections.abc").attr("Set");

    py_AdjacencyOuterDict.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(
        py_AdjacencyOuterDict.attr("__bases__"));
}

}