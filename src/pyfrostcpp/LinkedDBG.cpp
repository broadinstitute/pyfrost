#include "LinkedDBG.h"

namespace pyfrost {

void define_LinkedDBG(py::module& m) {
    py::class_<LinkedDBG<PyfrostCCDBG>>(m, "LinkedDBG")
        .def(py::init<PyfrostCCDBG*>())
        .def("add_links_from_seq", &LinkedDBG<PyfrostCCDBG>::addLinksFromSequence)
        .def("get_links", [] (LinkedDBG<PyfrostCCDBG>& self, Kmer const& kmer) {
            return self.getLinks(kmer);
        }, py::keep_alive<0, 1>())

        .def("__getitem__", [] (LinkedDBG<PyfrostCCDBG>& self, const Kmer& kmer) {
            return self.getLinks(kmer);
        }, py::keep_alive<0, 1>())
        .def("__len__", [] (LinkedDBG<PyfrostCCDBG>& self) {
            return self.numTrees();
        })
        .def("__iter__", [] (LinkedDBG<PyfrostCCDBG>& self) {
            return py::make_key_iterator(self.getJunctionTrees().begin(), self.getJunctionTrees().end());
        }, py::keep_alive<0, 1>());

}

}
