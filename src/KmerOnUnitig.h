#include <pybind11/pybind11.h>
#include "Pyfrost.h"

#ifndef PYFROST_KMERONUNITIG_H
#define PYFROST_KMERONUNITIG_H

namespace py = pybind11;

namespace pyfrost {


void define_KmerOnUnitig(py::module& m) {
    py::class_<PyfrostColoredUMap>(m, "KmerOnUnitig")
        .def_property_readonly("head", &PyfrostColoredUMap::getMappedHead, "The k-mer at the beginning of this unitig.")
        .def_property_readonly("tail", &PyfrostColoredUMap::getMappedTail, "The k-mer at the end of this unitig.")

        .def_readonly("pos", &PyfrostColoredUMap::dist, "Position of this k-mer on the unitig.")

        .def("__str__", [] (PyfrostColoredUMap const& self) {
            if(self.strand) {
                return self.referenceUnitigToString();
            } else {
                return reverse_complement(self.referenceUnitigToString());
            }
        })

        .def("__repr__", [] (PyfrostColoredUMap const& self) {
            stringstream repr;
            repr << "<KmerOnUnitig Unitig=";
            if(self.size > (2*Kmer::k)) {
                repr << self.getMappedHead().toString() << "..." << self.getMappedTail().toString();
            } else {
                if(self.strand) {
                    repr << self.referenceUnitigToString();
                } else {
                    repr << reverse_complement(self.referenceUnitigToString());
                }
            }
            repr << " Pos=" << self.dist << ">";

            return repr.str();
        })

        .def("__len__", [](PyfrostColoredUMap const& self) {
            return self.size;
        }, "Number of k-mers this unitig is composed of.")
        .def("__bool__", [] (PyfrostColoredUMap const& self) { return !self.isEmpty; });
}

}




#endif //PYFROST_KMERONUNITIG_H
