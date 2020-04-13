#ifndef PYFROST_PYFROST_H
#define PYFROST_PYFROST_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <ColoredCDBG.hpp>

#include "UnitigDataDict.h"

namespace py = pybind11;

namespace pyfrost {

using PyfrostCCDBG = ColoredCDBG<UnitigDataDict>;

enum class Strand : uint8_t {
    REVERSE,
    FORWARD
};

}

#endif //PYFROST_PYFROST_H
