#include <pybind11/pybind11.h>
#include "Pyfrost.h"
#include "UnitigDataProxy.h"

#ifndef PYFROST_KMERONUNITIG_H
#define PYFROST_KMERONUNITIG_H

namespace py = pybind11;

namespace pyfrost {

void define_UnitigMapping(py::module& m) {
    py::enum_<Strand>(m, "Strand")
        .value("REVERSE", Strand::REVERSE, "Reverse orientation")
        .value("FORWARD", Strand::FORWARD, "Forward orientation");

    py::class_<PyfrostColoredUMap>(m, "UnitigMapping")
        // Unitig head and tail k-mers
        .def_property_readonly("head", [] (PyfrostColoredUMap const& self) {
            return self.strand ? self.getUnitigHead() : self.getUnitigTail().twin();
            }, "The k-mer at the beginning of this unitig, in the same orientation as the mapped sequence.")
        .def_property_readonly("tail", [] (PyfrostColoredUMap const& self) {
            return self.strand ? self.getUnitigTail() : self.getUnitigHead().twin();
            }, "The k-mer at the end of this unitig, in the same orientation as the mapped sequence.")

        .def("get_full_mapping", &PyfrostColoredUMap::mappingToFullUnitig,
            "Get a new UnitigMapping object where the mapping represents the full unitig (instead of a single k-mer).")

        .def_property_readonly("data", [] (PyfrostColoredUMap const& self) {
            // Creating a new object every time `data` is accessed for now, maybe something smarter in the future
            return UnitigDataProxy(self);
        }, py::return_value_policy::move)

        .def("__str__", &PyfrostColoredUMap::mappedSequenceToString)

        .def("__repr__", [] (PyfrostColoredUMap const& self) {
            stringstream repr;

            if(self.isEmpty) {
                return std::string("<UnitigMapping EMPTY>");
            }

            repr << "<UnitigMapping Unitig=";
            if(self.size > (2*self.getGraph()->getK())) {
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

            if(self.len > self.getGraph()->getK()) {
                auto mapped_head = self.getMappedHead();
                auto mapped_tail = self.getMappedTail();

                repr << " MapStr=" << mapped_head.toString() << "..." << mapped_tail.toString();
            } else {
                repr << " MapStr=" << self.mappedSequenceToString();
            }

            repr << " MapPos=" << self.dist;
            repr << " MapLen=" << UnitigDataProxy(self).mappedSequenceLength();
            repr << " UnitigSize=" << self.size;
            repr << " Strand=" << (self.strand ? "forward" : "reverse") << ">";

            return repr.str();
        })

        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__hash__", [] (PyfrostColoredUMap const& self) {
            PyfrostColoredUMap::hasher hasher;

            return hasher(self);
        })

        .def("__bool__", [] (PyfrostColoredUMap const& self) { return !self.isEmpty; });
}

}




#endif //PYFROST_KMERONUNITIG_H
