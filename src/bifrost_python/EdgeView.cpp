#include "BifrostDiGraph.h"
#include "EdgeView.h"

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
        n2.strand ? Strand::FORWARD : Strand::REVERSE
    );

    return meta;
}

void define_EdgeView(py::module& m) {
    define_EdgeDataView<void>(m, "WholeGraph");
    define_EdgeDataView<py::iterable>(m, "PyCollection");

    auto py_OutEdgeView = py::class_<EdgeView<true>>(m, "OutEdgeView")
        .def("__getitem__", py::overload_cast<py::tuple const&>(&EdgeView<true>::get), py::is_operator())
        .def("__contains__", &EdgeView<true>::contains, py::is_operator())
        .def("__len__", &EdgeView<true>::size)

        .def_static("_from_iterable", [] (EdgeView<true>& self, py::iterable const& it) {
            return py::set(it);
        })

        .def("__iter__", [] (EdgeView<true> const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        // Make it callable and return a OutEdgeDataView if requested
        .def("__call__", [] (EdgeView<true> const& self) {
            // No arguments given, just iterate over all edges
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        .def("__call__", [] (EdgeView<true>& self, bool data, py::object const& default_value) {
            // All edges, but with data
            NodeIterable<void> iterable(self.getGraph());
            py::object py_data = py::bool_(data);
            return EdgeDataView<void, true>(self, iterable, py_data, default_value);
        }, py::keep_alive<0, 1>(), py::arg("data"), py::arg("default") = py::none())

        .def("__call__", [] (EdgeView<true>& self, std::string const& data, py::object const& default_value) {
            // All edges, but with data
            NodeIterable<void> iterable(self.getGraph());
            py::str py_data = py::str(data);
            return EdgeDataView<void, true>(self, iterable, py_data, default_value);
        }, py::keep_alive<0, 1>(), py::arg("data"), py::arg("default") = py::none())

        .def("__call__", [] (EdgeView<true>& self, py::iterable const& nbunch, py::object const& data,
                             py::object const& default_value) {
            // Now with potential data or given nodes
            NodeIterable<py::iterable> iterable(self.getGraph(), nbunch);
            return EdgeDataView<py::iterable, true>(self, iterable, data, default_value);
        }, py::arg("nbunch"), py::arg("data") = false, py::arg("default") = py::none(),
            py::keep_alive<0, 1>(), py::keep_alive<0, 2>());

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    auto Set = py::module::import("collections.abc").attr("Set");
    py_OutEdgeView.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(py_OutEdgeView.attr("__bases__"));

    auto py_InEdgeView = py::class_<EdgeView<false>>(m, "InEdgeView")
        .def("__getitem__", py::overload_cast<py::tuple const&>(&EdgeView<false>::get), py::is_operator())
        .def("__contains__", &EdgeView<false>::contains, py::is_operator())
        .def("__len__", &EdgeView<false>::size)

        .def_static("_from_iterable", [] (EdgeView<true>& self, py::iterable const& it) {
            return py::set(it);
        })

        .def("__iter__", [] (EdgeView<false> const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        // Make it callable and return an InEdgeDataView if requested
        .def("__call__", [] (EdgeView<false> const& self) {
            // No arguments given, just iterate over all edges
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        .def("__call__", [] (EdgeView<false>& self, bool data, py::object const& default_value) {
            // All edges, but with data
            NodeIterable<void> iterable(self.getGraph());
            py::object py_data = py::bool_(data);
            return EdgeDataView<void, false>(self, iterable, py_data, default_value);
        }, py::keep_alive<0, 1>(), py::arg("data"), py::arg("default") = py::none())

        .def("__call__", [] (EdgeView<false>& self, std::string const& data, py::object const& default_value) {
            // All edges, but with data
            NodeIterable<void> iterable(self.getGraph());
            py::str py_data = py::str(data);
            return EdgeDataView<void, false>(self, iterable, py_data, default_value);
        }, py::keep_alive<0, 1>(), py::arg("data"), py::arg("default") = py::none())

        .def("__call__", [] (EdgeView<false>& self, py::iterable const& nbunch, py::object const& data,
                             py::object const& default_value) {
             // Now with potential data or given nodes
             NodeIterable<py::iterable> iterable(self.getGraph(), nbunch);
             return EdgeDataView<py::iterable, false>(self, iterable, data, default_value);
         }, py::arg("nbunch"), py::arg("data") = false, py::arg("default") = py::none(),
         py::keep_alive<0, 1>(), py::keep_alive<0, 2>());

    py_InEdgeView.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(py_InEdgeView.attr("__bases__"));
}

}
