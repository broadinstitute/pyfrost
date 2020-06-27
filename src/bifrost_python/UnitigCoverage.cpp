#include "UnitigCoverage.h"

namespace pyfrost {

void define_UnitigCoverageProxy(py::module& m)
{
    py::class_<UnitigCoverageProxy>(m, "UnitigCoverageProxy")
        .def("__getitem__", &UnitigCoverageProxy::getCoverage, py::is_operator())
        .def("__len__", &UnitigCoverageProxy::size, py::is_operator())

        .def("__iter__", [](UnitigCoverageProxy& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());
}

}
