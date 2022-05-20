#include "Kmer.h"

namespace pyfrost {


Kmer to_kmer(PyfrostColoredUMap const& obj, bool rev_compl) {
    if(rev_compl) {
        return obj.getMappedTail().twin();
    } else {
        return obj.getMappedHead();
    }
}

Kmer to_kmer(py::handle const& obj, bool rev_compl) {
    if(obj) {
        if(py::isinstance<Kmer>(obj)) {
            return py::cast<Kmer>(obj);
        } else if (py::isinstance<py::str>(obj)) {
            std::string kmer = py::cast<std::string>(obj);
            if(kmer.size() == Kmer::k) {
                auto kmer_obj = Kmer(kmer.c_str());
                if(rev_compl) {
                    return kmer_obj.twin();
                } else {
                    return kmer_obj;
                }
            }
        }
    }

    Kmer kmer;
    kmer.set_empty();

    return kmer;
}

Kmer to_kmer(py::object const& obj, bool rev_compl) {
    if(py::isinstance<Kmer>(obj)) {
        return py::cast<Kmer>(obj);
    } else if (py::isinstance<py::str>(obj)) {
        std::string kmer = py::cast<std::string>(obj);
        if(kmer.size() == Kmer::k) {
            auto kmer_obj = Kmer(kmer.c_str());
            if(rev_compl) {
                return kmer_obj.twin();
            } else {
                return kmer_obj;
            }
        }
    }

    Kmer kmer;
    kmer.set_empty();

    return kmer;
}

Kmer to_kmer(Kmer const& kmer, bool rev_compl) {
    if(rev_compl) {
        return kmer.twin();
    } else {
        return kmer;
    }
}

Kmer to_kmer(char const* kmer, bool rev_compl) {
    if(strlen(kmer) == Kmer::k) {
        auto kmer_obj = Kmer(kmer);
        if(rev_compl) {
            return kmer_obj.twin();
        } else {
            return kmer_obj;
        }
    }

    Kmer kmer_obj;
    kmer_obj.set_empty();

    return kmer_obj;
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
            if(index < 0) {
                index = Kmer::k + index;
            }

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

        .def("to_binary", [] (Kmer const& self) {
            stringstream stream;
            self.write(stream);

            std::string kmer_binary = stream.str();
            return py::bytes(kmer_binary);
        })

        // Other functions
        .def("__str__", [] (Kmer const& self) { return self.toString(); })
        .def("__repr__", [] (Kmer const& self) { return "<Kmer '" + self.toString() + "'>"; })
        .def("__bool__", [] (Kmer const& self) {
            return !is_kmer_empty(self);
        })
        .def("__hash__", [] (Kmer const& self) { return self.hash(); })
        .def("hash", &Kmer::hash, "Get the hash value for this k-mer", py::arg("seed") = 0);

    /**
     * Class masking as a function, can't use a single function here
     * because we need to keep the string to k-merize in memory while
     * iterating over it.
     */
    py::class_<Kmerizer>(m, "kmerize_str")
        .def(py::init<std::string const&>())
        .def("__iter__", [] (Kmerizer const& self) {
            return py::make_key_iterator<py::return_value_policy::copy>(self.begin(), self.end());
        }, py::keep_alive<0, 1>());


    // two bits are reserved for k-mer metadata
    m.attr("max_k") = MAX_KMER_SIZE - 1;

    m.def("set_k", &Kmer::set_k,
          "Set the k-mer size for all `Kmer` objects. Only use at the beginning of the program!");
}

}
