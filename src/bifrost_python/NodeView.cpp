#include "NodeView.h"

namespace pyfrost {

void define_NodeView(py::module& m) {
    auto py_NodeDataView = py::class_<NodeDataView>(m, "NodeDataView")
        .def("__getitem__", py::overload_cast<char const*>(&NodeDataView::get), py::is_operator())
        .def("__getitem__", py::overload_cast<Kmer const&>(&NodeDataView::get), py::is_operator())

        .def("__contains__", py::overload_cast<char const*>(&NodeDataView::contains), py::is_operator())
        .def("__contains__", py::overload_cast<Kmer const&>(&NodeDataView::contains), py::is_operator())

        .def("__len__", &NodeDataView::numNodes, py::is_operator())

        .def("__iter__", [] (NodeDataView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

            // Required for Set mixin (see collections.abc Python documentation)
        .def_static("_from_iterable", [] (py::iterable const& iterator) {
            return py::set(iterator);
        });

    // Hack to let our NodeDataView inherit from the collections.abc.Set Python mixin
    auto Set = py::module::import("collections.abc").attr("Set");
    py_NodeDataView.attr("__bases__") = py::make_tuple(Set).attr("__add__")(py_NodeDataView.attr("__bases__"));

    auto py_NodeView = py::class_<NodeView>(m, "NodeView")
        // Access nodes with [] operator overloading
        .def("__getitem__", py::overload_cast<char const*>(&NodeView::findNode), py::is_operator())
        .def("__getitem__", py::overload_cast<Kmer const&>(&NodeView::findNode), py::is_operator())

        .def("__contains__", py::overload_cast<char const*>(&NodeView::contains), py::is_operator())
        .def("__contains__", py::overload_cast<Kmer const&>(&NodeView::contains), py::is_operator())

        .def("__len__", &NodeView::numNodes, py::is_operator())

            // Call without parameters just returns the node iterator
        .def("__call__", [] (NodeView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

            // With data arguments we return a NodeDataView
        .def("__call__", [] (NodeView& self, py::object const& data, py::object const& default_value) {
                 try {
                     bool data_bool = py::cast<bool>(data);

                     return new NodeDataView(self, data_bool);
                 } catch(py::cast_error const&) {
                     return new NodeDataView(self, true, data, default_value);
                 }
             }, py::arg("data") = false, py::arg("default") = py::none(),
             py::keep_alive<0, 1>())

            // Iterate over all nodes
        .def("__iter__", [](NodeView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

            // Required for Set mixin (see collections.abc Python documentation)
        .def_static("_from_iterable", [] (py::iterable const& iterator) {
            return py::set(iterator);
        });

    // NodeView inherits from both Mapping and Set
    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    py_NodeView.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(py_NodeView.attr("__bases__"));
}

}
