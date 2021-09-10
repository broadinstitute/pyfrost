#include "LinkDB.h"

namespace pyfrost {

void define_LinkDB(py::module& m) {
    // Simple abstract class
    py::class_<LinkDB>(m, "LinkDB");
}

}
