#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <Kmer.hpp>

namespace py = pybind11;

#ifndef PYFROST_KMER_H
#define PYFROST_KMER_H

namespace pyfrost {

void define_Kmer(py::module &m) {
    py::class_<Kmer>(m, "Kmer")
        // Constructors
        .def(py::init<>())
        .def(py::init<Kmer const&>())
        .def(py::init<char const*>())

        // Operator overloading
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self < py::self)
        .def("__getitem__", [] (const Kmer& self, int index) {
                if(index >= Kmer::k || index < 0) {
                    throw std::out_of_range("Index is out of range");
                }

                return self.getChar(index);
            }, py::is_operator(), "Get base at position")

        // DNA sequence operations
        .def("forward_base", &Kmer::forwardBase, "Shift the k-mer to the left with one base and append the given base.")
        .def("backward_base", &Kmer::forwardBase, "Shift the k-mer to the right with one base and prepend the given base.")
        .def("twin", &Kmer::twin, "Get the reverse complement of this k-mer.")
        .def("rep", &Kmer::rep, "Get the canonical k-mer.")

        // Other functions
        .def("__str__", [] (const Kmer& self) { return self.toString(); })
        .def("__repr__", [] (const Kmer& self) { return "<Kmer '" + self.toString() + "'>"; })
        .def("__hash__", [] (const Kmer& self) { return self.hash(); })
        .def("hash", &Kmer::hash, "Get the hash value for this k-mer", py::arg("seed") = 0);


    m.attr("max_k") = MAX_KMER_SIZE;

    m.def("set_k", &Kmer::set_k, "Set the k-mer size for all `Kmer` objects. Only use at the beginning of the program!");
}

}

#endif //PYFROST_KMER_H
