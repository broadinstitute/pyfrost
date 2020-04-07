#include <ColoredCDBG.hpp>

#include "UnitigDataDict.h"

#ifndef PYFROST_PYFROST_H
#define PYFROST_PYFROST_H

namespace pyfrost {

using PyfrostCCDBG = ColoredCDBG<UnitigDataDict>;

enum class Strand : uint8_t {
    REVERSE,
    FORWARD
};

}

#endif //PYFROST_PYFROST_H
