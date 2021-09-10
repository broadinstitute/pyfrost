#include "LinkAnnotator.h"


namespace pyfrost {

void define_LinkAnnotator(py::module& m) {
    py::class_<LinkAnnotator<PyfrostCCDBG>>(m, "LinkAnnotator")
        .def(py::init())
        .def("add_links_from_sequence",
             py::overload_cast<PyfrostCCDBG&, LinkDB&, string const&, bool>(&LinkAnnotator<PyfrostCCDBG>::addLinksFromSequence),
             py::arg("graph"), py::arg("db"), py::arg("sequence"), py::arg("keep_nodes") = false);

    py::class_<RefLinkAnnotator<PyfrostCCDBG>, LinkAnnotator<PyfrostCCDBG>>(m, "RefLinkAnnotator")
        .def(py::init());

    m.def("add_links_from_fasta", &addLinksFromFasta<PyfrostCCDBG>);
}

void define_MappingResult(py::module& m) {
    py::class_<MappingResult>(m, "MappingResult")
        .def_readonly("start_unitig", &MappingResult::start_unitig)
        .def_readonly("end_unitig", &MappingResult::end_unitig)
        .def_readonly("mapping_start", &MappingResult::mapping_start)
        .def_readonly("mapping_end", &MappingResult::mapping_end)
        .def_readonly("num_junctions", &MappingResult::num_junctions)
        .def_readonly("unitig_visits", &MappingResult::unitig_visits)
        .def("matching_kmers", [] (MappingResult& self) {
            // Make copy and return as NumPy array
            vector<uint8_t> to_return(self.matches);

            return as_pyarray<vector<uint8_t>>(move(to_return));
        })
        .def("__str__", [] (MappingResult& self) {
            stringstream sstream;
            sstream << self.start_unitig.toString() << '\t'
                    << self.end_unitig.toString() << '\t'
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
