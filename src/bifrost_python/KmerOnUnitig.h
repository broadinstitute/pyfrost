#include <pybind11/pybind11.h>
#include "Pyfrost.h"

#ifndef PYFROST_KMERONUNITIG_H
#define PYFROST_KMERONUNITIG_H

namespace py = pybind11;

namespace pyfrost {


void define_KmerOnUnitig(py::module& m) {
    py::class_<PyfrostColoredUMap>(m, "KmerOnUnitig")
        .def_property_readonly("head", [] (PyfrostColoredUMap const& self) {
            return self.strand ? self.getUnitigHead() : self.getUnitigTail().twin();
            }, "The k-mer at the beginning of this unitig.")
        .def_property_readonly("tail", [] (PyfrostColoredUMap const& self) {
            return self.strand ? self.getUnitigTail() : self.getUnitigHead().twin();
            }, "The k-mer at the end of this unitig.")

        .def_readonly("pos", &PyfrostColoredUMap::dist, "Position of this k-mer on the unitig.")

        .def("__str__", [] (PyfrostColoredUMap const& self) {
            if(self.strand) {
                return self.referenceUnitigToString();
            } else {
                return reverse_complement(self.referenceUnitigToString());
            }
        })

        .def_readonly("__len__", &PyfrostColoredUMap::size, "Length of the unitig sequence.")

        .def("__repr__", [] (PyfrostColoredUMap const& self) {
            stringstream repr;

            repr << "<KmerOnUnitig Unitig=";
            if(self.size > (2*Kmer::k)) {
                auto head = self.strand ? self.getUnitigHead() : self.getUnitigTail().twin();
                auto tail = self.strand ? self.getUnitigTail() : self.getUnitigHead().twin();

                repr << head.toString() << "..." << tail.toString();
            } else {
                auto unitig_str = self.referenceUnitigToString();
                if(!self.strand) {
                    unitig_str = reverse_complement(unitig_str);
                }

                repr << unitig_str;
            }
            repr << " Pos=" << self.dist;
            repr << " " << "Strand=" << (self.strand ? "forward" : "reverse") << ">";

            return repr.str();
        })

        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__hash__", [] (PyfrostColoredUMap const& self) { return self.getUnitigHead().hash(); })

        .def("__bool__", [] (PyfrostColoredUMap const& self) { return !self.isEmpty; });
}

}




#endif //PYFROST_KMERONUNITIG_H
