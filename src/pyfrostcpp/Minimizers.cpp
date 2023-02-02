#include "Minimizers.h"
#include "Kmer.hpp"


namespace pyfrost {

void define_Minimizer(py::module& m) {
    py::class_<Minimizer>(m, "Minimizer")
        // Constructors
        .def(py::init<>())
        .def(py::init<Minimizer const&>())
        .def(py::init<char const*>())
        .def_static("from_bytes", &Minimizer::fromBytes, "Create minimizer from bytes.")

        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self < py::self)

        .def("twin", &Minimizer::twin, "Get the reverse complement of this Minimizer")
        .def("rep", &Minimizer::rep, "Get the canonical minimizer")

        .def("__bytes__", [] (Minimizer const& self) {
                return py::bytes(self.asBytes());
        }, "Get binary representation")
        .def("__str__", py::overload_cast<>(&Minimizer::toString, py::const_), "Convert to string")
        .def("__repr__", [] (Minimizer const& self) {
            stringstream ss;
            ss << "<Minimizer '" << self.toString() << "'>";

            return ss.str();
        })
        .def("__len__", [] (Minimizer const& self) { return Minimizer::g; })
        .def("__bool__", [] (Minimizer const& self) {
            return !(self.isEmpty() || self.isDeleted());
        })
        .def("__hash__", [] (Minimizer const& self) {
            RepHash hasher(Minimizer::g);
            hasher.init(self.toString().c_str());

            return hasher.hash();
        });
}

void define_MinHashIterator(py::module& m) {
    py::class_<MinHashIterator>(m, "minhash_iter")
        .def(py::init<std::string const&>())
        .def(py::init<std::string const&, size_t>())
        .def(py::init<std::string const&, size_t, size_t>())

        .def("__iter__", [] (MinHashIterator const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.begin(), self.end());
        }, py::keep_alive<0, 1>());
}

void define_MinHashResult(py::module& m) {
    py::class_<minHashResult>(m, "MinHashResult")
        .def_readonly("hash", &minHashResult::hash)
        .def_readonly("pos", &minHashResult::pos)

        .def("__repr__", [] (minHashResult const& self) {
            stringstream ss;

            ss << "<MinHashResult hash=" << self.hash << " pos=" << self.pos << ">";
            return ss.str();
        });
}


}
