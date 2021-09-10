#ifndef PYFROST_PYFROST_H
#define PYFROST_PYFROST_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include <ColoredCDBG.hpp>

#include "UnitigDataDict.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<std::string>);
PYBIND11_MAKE_OPAQUE(std::unordered_map<Kmer, size_t>);


namespace pyfrost {

/**
 * Utility function to expose any C++ sequence container as numpy array without copying.
 *
 * Source: https://github.com/pybind/pybind11/issues/1042
 */
template <typename Sequence,
          typename = std::enable_if_t<std::is_rvalue_reference<Sequence&&>::value>>
inline py::array_t<typename Sequence::value_type> as_pyarray(Sequence&& seq) {
    auto size = seq.size();
    auto data = seq.data();
    std::unique_ptr<Sequence> seq_ptr = std::make_unique<Sequence>(std::forward<Sequence>(seq));
    auto capsule = py::capsule(seq_ptr.get(), [](void *p) { std::unique_ptr<Sequence>(reinterpret_cast<Sequence*>(p)); });
    seq_ptr.release();
    return py::array(size, data, capsule);
}

using PyfrostCCDBG = ColoredCDBG<UnitigDataDict>;

enum class Strand : uint8_t {
    REVERSE,
    FORWARD
};

}

#endif //PYFROST_PYFROST_H
