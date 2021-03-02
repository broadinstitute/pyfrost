#ifndef PYFROST_KMER_H
#define PYFROST_KMER_H

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "UnitigDataDict.h"

#include <Kmer.hpp>

// Specialize hashing for Kmer objects (to be used with robin_hood unordered_map)
namespace std {

template<>
class hash<Kmer> {
public:
    size_t operator()(Kmer const& kmer) const {
        return kmer.hash();
    }
};

}

namespace py = pybind11;

namespace pyfrost {

/**
 * Convert any object to a K-mer. By default simply outputs an empty k-mer, but specializations exist to convert to
 * more meaningful values.
 *
 * This function is for example used by NodeIterator to convert any value in a given Python list to a potential node
 * k-mer.
 */
template <typename T>
Kmer to_kmer(T, bool rev_compl=false) {
    Kmer kmer;
    kmer.set_empty();

    return kmer;
}

Kmer to_kmer(PyfrostColoredUMap const& obj, bool rev_compl=false);

Kmer to_kmer(py::handle const& obj, bool rev_compl=false);

Kmer to_kmer(py::object const& obj, bool rev_compl=false);

Kmer to_kmer(Kmer const& kmer, bool rev_compl=false);

Kmer to_kmer(char const* kmer, bool rev_compl=false);

bool is_kmer_empty(Kmer const& kmer);

void define_Kmer(py::module &m);


class Kmerizer {
public:
    explicit Kmerizer(std::string const& _str) : str(_str) { }
    Kmerizer(Kmerizer const& o) = default;
    Kmerizer(Kmerizer&& o) = default;

    KmerIterator begin() const {
        return KmerIterator(str.c_str());
    }

    KmerIterator end() const {
        return KmerIterator();
    }

private:
    std::string str;
};

}

#endif //PYFROST_KMER_H
