#include "AdjacencyProxy.h"

namespace pyfrost {

void define_AdjacencyProxy(py::module& m) {
    py::enum_<AdjacencyType>(m, "AdjacencyType")
        .value("SUCCESSORS", AdjacencyType::SUCCESSORS)
        .value("PREDECESSORS", AdjacencyType::PREDECESSORS);

    py::class_<AdjacencyViewBase>(m, "AdjacencyViewBase")
        .def("__len__", &AdjacencyViewBase::numNeigbors, py::is_operator())

        .def("__contains__", py::overload_cast<char const*>(&AdjacencyViewBase::contains), py::is_operator())
        .def("__contains__", py::overload_cast<Kmer const&>(&AdjacencyViewBase::contains), py::is_operator())
        .def("__contains__", py::overload_cast<PyfrostColoredUMap const&>(&AdjacencyViewBase::contains), py::is_operator())

        .def("__getitem__", &AdjacencyViewBase::getEdgeDict);

    py::class_<SuccessorView, AdjacencyViewBase>(m, "SuccessorView")
        .def(py::init<PyfrostColoredUMap const&>())
        .def("__iter__", [] (SuccessorView const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    py::class_<PredecessorView, AdjacencyViewBase>(m, "PredecessorView")
        .def(py::init<PyfrostColoredUMap const&>())
        .def("__iter__", [] (PredecessorView const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    py::class_<AdjacencyProxy>(m, "AdjacencyProxy")
        .def("__len__", &AdjacencyProxy::numNodes, py::is_operator())

        .def("__contains__", py::overload_cast<PyfrostColoredUMap const&>(&AdjacencyProxy::contains))
        .def("__contains__", py::overload_cast<Kmer const&>(&AdjacencyProxy::contains))
        .def("__contains__", py::overload_cast<char const*>(&AdjacencyProxy::contains))

        .def("__getitem__", py::overload_cast<PyfrostColoredUMap const&>(&AdjacencyProxy::getView), py::is_operator())
        .def("__getitem__", py::overload_cast<Kmer const&>(&AdjacencyProxy::getView), py::is_operator())
        .def("__getitem__", py::overload_cast<char const*>(&AdjacencyProxy::getView), py::is_operator())

        .def("__iter__", [] (AdjacencyProxy const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.begin(), self.end());
        }, py::keep_alive<0, 1>());
}

}