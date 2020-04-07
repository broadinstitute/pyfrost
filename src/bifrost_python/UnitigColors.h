#include <pybind11/pybind11.h>
#include <ColorSet.hpp>

namespace py = pybind11;

#ifndef PYFROST_UNITIGCOLORS_H
#define PYFROST_UNITIGCOLORS_H

namespace pyfrost {

void define_UnitigColors(py::module& m) {
    py::class_<UnitigColors>(m, "UnitigColors");
        /*.def("__iter__", [] (UnitigColors const& self) {
            return py::make_iterator(self.begin(), self.end());
        })*/
}

}



#endif //PYFROST_UNITIGCOLORS_H
