#include "LinkAnnotator.h"


namespace pyfrost {

void define_LinkAnnotator(py::module& m) {
    py::class_<LinkAnnotator<PyfrostCCDBG>>(m, "LinkAnnotator")
        .def(py::init<PyfrostCCDBG*, LinkDB*>())
        .def_property("max_link_length",
                      &LinkAnnotator<PyfrostCCDBG>::getMaxLinkLength, &LinkAnnotator<PyfrostCCDBG>::setMaxLinkLength)
        .def("add_links_from_sequence",
             py::overload_cast<string const&, bool>(&LinkAnnotator<PyfrostCCDBG>::addLinksFromSequence),
             py::arg("sequence"), py::arg("keep_nodes") = false)
         .def("add_links_from_path", &LinkAnnotator<PyfrostCCDBG>::addLinksFromPath);

    py::class_<ColorAssociatedAnnotator<PyfrostCCDBG>, LinkAnnotator<PyfrostCCDBG>>(m, "ColorAssociatedAnnotator")
        .def(py::init<PyfrostCCDBG*, LinkDB*>());

    m.def("add_links_from_fasta", &addLinksFromFasta<PyfrostCCDBG>, py::call_guard<py::gil_scoped_release>());
}

void define_MappingResult(py::module& m) {
    py::class_<MappingResult>(m, "MappingResult")
        .def_readonly("paths", &MappingResult::paths)
        .def("start_unitig", &MappingResult::start_unitig)
        .def("end_unitig", &MappingResult::end_unitig)
        .def_readonly("mapping_start", &MappingResult::mapping_start)
        .def_readonly("mapping_end", &MappingResult::mapping_end)
        .def_readonly("unitig_visits", &MappingResult::unitig_visits)
        .def("matching_kmers", [] (MappingResult& self) {
            // Make copy and return as NumPy array
            vector<uint8_t> to_return(self.matches);

            return as_pyarray<vector<uint8_t>>(move(to_return));
        })
        .def("__str__", [] (MappingResult& self) {
            stringstream sstream;
            sstream << self.start_unitig().toString() << '\t'
                    << self.end_unitig().toString() << '\t'
                    << self.matches.size() << '\t'
                    << self.mapping_start << '\t'
                    << self.mapping_end << '\t';

            size_t num_correct = 0;
            for(unsigned char match : self.matches) {
                if(match > 0) {
                    ++num_correct;
                }

                sstream << (match > 0 ? '1' : '0');
            }
            sstream.precision(2);
            sstream << '\t' << ((num_correct * 100.0) / self.matches.size());

            return sstream.str();
        });
}

}
