#include "AdjacencyProxy.h"

namespace pyfrost {

void define_AdjacencyProxy(py::module& m) {
    py::enum_<AdjacencyType>(m, "AdjacencyType")
        .value("SUCCESSORS", AdjacencyType::SUCCESSORS)
        .value("PREDECESSORS", AdjacencyType::PREDECESSORS);

    auto py_AdjacencyViewBase = py::class_<AdjacencyViewBase>(m, "AdjacencyViewBase")
        .def("__iter__", [] (AdjacencyViewBase const&) {
            throw py::type_error("Calling __iter__ on the abstract base class!");
        })
        .def("__len__", &AdjacencyViewBase::numNeigbors, py::is_operator())

        .def("__contains__", py::overload_cast<char const*>(&AdjacencyViewBase::contains), py::is_operator())
        .def("__contains__", py::overload_cast<Kmer const&>(&AdjacencyViewBase::contains), py::is_operator())

        .def("__getitem__", [] (AdjacencyViewBase const& self, const Kmer&) {
            return py::dict(); // TODO
        });

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    auto Set = py::module::import("collections.abc").attr("Set");
    py_AdjacencyViewBase.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(
        py_AdjacencyViewBase.attr("__bases__"));

    py::class_<SuccessorView, AdjacencyViewBase>(m, "SuccessorView")
        .def("__iter__", [] (SuccessorView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    py::class_<PredecessorView, AdjacencyViewBase>(m, "PredecessorView")
        .def("__iter__", [] (PredecessorView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    auto py_AdjacencyProxy = py::class_<AdjacencyProxy>(m, "AdjacencyProxy")
        .def("__len__", &AdjacencyProxy::numNodes, py::is_operator())

        .def("__contains__", py::overload_cast<Kmer const&>(&AdjacencyProxy::contains))
        .def("__contains__", py::overload_cast<char const*>(&AdjacencyProxy::contains))

        .def("__getitem__", py::overload_cast<Kmer const&>(&AdjacencyProxy::getView), py::is_operator())
        .def("__getitem__", py::overload_cast<char const*>(&AdjacencyProxy::getView), py::is_operator())

        .def("__iter__", [] (AdjacencyProxy const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    py_AdjacencyProxy.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(
        py_AdjacencyProxy.attr("__bases__"));
}

}