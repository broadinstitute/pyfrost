#include "UnitigColors.h"

namespace pyfrost {

void define_UnitigColors(py::module &m) {
    py::class_<UnitigColorsProxy>(m, "UnitigColors")
        .def("__contains__", &UnitigColorsProxy::contains)
        .def("__len__", &UnitigColorsProxy::size)
        .def("__iter__", [](UnitigColorsProxy const &self) {
            return py::make_iterator(self.begin(), self.end());
        })

        .def("num_kmers_with_color", &UnitigColorsProxy::numKmersWithColor);
}

}
