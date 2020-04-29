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

void define_EdgeView(py::module& m) {
    py::class_<EdgeDataView<true>>(m, "OutEdgeDataView")
        .def("__getitem__", &EdgeDataView<true>::get, py::is_operator())
        .def("__contains__", &EdgeDataView<true>::contains, py::is_operator())
        .def("__len__", &EdgeDataView<true>::size, py::is_operator())
        .def("__iter__", [] (EdgeDataView<true> const& self) {
            if(!self.hasNodeBunch()) {
                return py::make_iterator(self.all_begin(), self.all_end());
            } else {
                return py::make_iterator(self.nbunch_begin(), self.nbunch_end());
            }
        }, py::keep_alive<0, 1>());

    py::class_<EdgeDataView<false>>(m, "InEdgeDataView")
        .def("__getitem__", &EdgeDataView<false>::get, py::is_operator())
        .def("__contains__", &EdgeDataView<false>::contains, py::is_operator())
        .def("__len__", &EdgeDataView<false>::size, py::is_operator())
        .def("__iter__", [] (EdgeDataView<false> const& self) {
            if(!self.hasNodeBunch()) {
                return py::make_iterator(self.all_begin(), self.all_end());
            } else {
                return py::make_iterator(self.nbunch_begin(), self.nbunch_end());
            }
        }, py::keep_alive<0, 1>());

    auto py_OutEdgeView = py::class_<EdgeView<true>>(m, "OutEdgeView")
        .def("__getitem__", py::overload_cast<py::tuple const&>(&EdgeView<true>::get), py::is_operator())
        .def("__contains__", &EdgeView<true>::contains, py::is_operator())
        .def("__len__", &EdgeView<true>::size)
        .def("__iter__", [] (EdgeView<true> const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        // Make it callable and return a OutEdgeDataView if requested
        .def("__call__", [] (EdgeView<true> const& self) {
            // No arguments given, just iterate over edges
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())
        .def("__call__", [] (EdgeView<true>& self, py::iterable const& nbunch, py::object const& data,
                             py::object const& default_value) {
            // Now with potential data or given nodes
            return new EdgeDataView<true>(self, nbunch, data, default_value);
        }, py::arg("nbunch") = py::none(), py::arg("data") = false, py::arg("default") = py::none(),
            py::keep_alive<0, 1>())
        .def("__call__", [] (EdgeView<true>& self, py::kwargs const& kwargs) {
            auto nbunch = kwargs.contains("nbunch") ? kwargs["nbunch"].cast<py::object>() : py::none();
            auto data = kwargs.contains("data") ? kwargs["data"].cast<py::object>() : py::cast<bool>(false);
            auto default_value = kwargs.contains("default") ? kwargs["default"].cast<py::object>() : py::none();

            return new EdgeDataView<true>(self, nbunch, data, default_value);
        }, py::keep_alive<0, 1>());

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    auto Set = py::module::import("collections.abc").attr("Set");
    py_OutEdgeView.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(py_OutEdgeView.attr("__bases__"));

    auto py_InEdgeView = py::class_<EdgeView<false>>(m, "InEdgeView")
        .def("__getitem__", py::overload_cast<py::tuple const&>(&EdgeView<false>::get), py::is_operator())
        .def("__contains__", &EdgeView<false>::contains, py::is_operator())
        .def("__len__", &EdgeView<false>::size)
        .def("__iter__", [] (EdgeView<false> const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        // Make it callable and return an InEdgeDataView if requested
        .def("__call__", [] (EdgeView<false> const& self) {
            // No arguments given, just iterate over edges
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())
        .def("__call__", [] (EdgeView<false>& self, py::iterable const& nbunch, py::object const& data,
                             py::object const& default_value) {
             // Now with potential data or given nodes
             return new EdgeDataView<false>(self, nbunch, data, default_value);
         }, py::arg("nbunch") = py::none(), py::arg("data") = false, py::arg("default") = py::none(),
         py::keep_alive<0, 1>())
        .def("__call__", [] (EdgeView<false>& self, py::kwargs const& kwargs) {
            auto nbunch = kwargs.contains("nbunch") ? kwargs["nbunch"].cast<py::object>() : py::none();
            auto data = kwargs.contains("data") ? kwargs["data"].cast<py::object>() : py::cast<bool>(false);
            auto default_value = kwargs.contains("default") ? kwargs["default"].cast<py::object>() : py::none();

            return new EdgeDataView<false>(self, nbunch, data, default_value);
        }, py::keep_alive<0, 1>());

    py_InEdgeView.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(py_InEdgeView.attr("__bases__"));

}

}
