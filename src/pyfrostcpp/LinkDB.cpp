#include "LinkDB.h"

namespace pyfrost {

void define_LinkDB(py::module& m) {
    // Simple abstract class
    py::class_<LinkDB>(m, "LinkDB");

    m.def("add_links_from_sequence", &addLinksFromSequence<PyfrostCCDBG>);
    m.def("add_links_from_files", &addLinksFromFile<PyfrostCCDBG>);
}

}
