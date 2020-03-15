#include <pybind11/pybind11.h>
#include "Pyfrost.h"

#ifndef PYFROST_KMERONUNITIG_H
#define PYFROST_KMERONUNITIG_H

namespace py = pybind11;

namespace pyfrost {


void define_KmerOnUnitig(py::module& m) {
    py::class_<PyfrostColoredUMap>(m, "KmerOnUnitig")
        .def_property_readonly("head", &PyfrostColoredUMap::getUnitigHead, "The k-mer at the beginning of this unitig.")
        .def_property_readonly("tail", &PyfrostColoredUMap::getUnitigTail, "The k-mer at the end of this unitig.")

        .def("__str__", &PyfrostColoredUMap::referenceUnitigToString)
        .def("__bool__", [] (PyfrostColoredUMap const& self) { return !self.isEmpty; });

}

}




#endif //PYFROST_KMERONUNITIG_H
