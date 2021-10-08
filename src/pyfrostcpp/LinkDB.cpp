#include "LinkDB.h"

namespace pyfrost {

void define_LinkDB(py::module& m) {
    // Simple abstract class
    py::class_<LinkDB>(m, "LinkDB")
        .def_property("color",[] (LinkDB& self) -> py::object {
                          return self.getColor() ? py::cast(*self.getColor()) : py::none();
                      }, &LinkDB::setColor);
}

}
