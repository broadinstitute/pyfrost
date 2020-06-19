#include "Kmer.h"

namespace pyfrost {


Kmer to_kmer(PyfrostColoredUMap const& obj) {
    return obj.getMappedHead();
}

Kmer to_kmer(py::handle const& obj) {
    if(obj) {
        if(py::isinstance<Kmer>(obj)) {
            return py::cast<Kmer>(obj);
        } else if (py::isinstance<py::str>(obj)) {
            std::string kmer = py::cast<std::string>(obj);
            if(kmer.size() == Kmer::k) {
                return Kmer(kmer.c_str());
            }
        }
    }

    Kmer kmer;
    kmer.set_empty();

    return kmer;
}

Kmer to_kmer(py::object const& obj) {
    if(py::isinstance<Kmer>(obj)) {
        return py::cast<Kmer>(obj);
    } else if (py::isinstance<py::str>(obj)) {
        std::string kmer = py::cast<std::string>(obj);
        if(kmer.size() == Kmer::k) {
            return Kmer(kmer.c_str());
        }
    }

    Kmer kmer;
    kmer.set_empty();

    return kmer;
}

Kmer to_kmer(Kmer const& kmer) {
    return kmer;
}

Kmer to_kmer(char const* kmer) {
    return Kmer(kmer);
}

bool is_kmer_empty(Kmer const& kmer) {
    Kmer empty;
    empty.set_empty();

    return kmer == empty;
}

/**
 * Define the Kmer class in Python
 *
 * While the original Bifrost K-mer is a mutable object, our Python counterpart is inmutable to be used for keys in
 * dicts, sets etc.
 */
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
        .def("__bool__", [] (const Kmer& self) {
            Kmer empty;
            empty.set_empty();

            return self != empty;
        })
        .def("__hash__", [] (const Kmer& self) { return self.hash(); })
        .def("hash", &Kmer::hash, "Get the hash value for this k-mer", py::arg("seed") = 0);


    // two bits are reserved for k-mer metadata
    m.attr("max_k") = MAX_KMER_SIZE - 1;

    m.def("set_k", &Kmer::set_k, "Set the k-mer size for all `Kmer` objects. Only use at the beginning of the program!");
}

}
