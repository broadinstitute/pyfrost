#ifndef PYFROST_UNITIGCOLORS_H
#define PYFROST_UNITIGCOLORS_H

#include "Pyfrost.h"

namespace py = pybind11;

namespace pyfrost {

class UnitigColorsProxy {
public:
    explicit UnitigColorsProxy(PyfrostColoredUMap const& unitig) : unitig(unitig) { }

    UnitigColorsProxy(UnitigColorsProxy const& o) = default;
    UnitigColorsProxy(UnitigColorsProxy&& o) = default;

    UnitigColorsProxy& operator=(UnitigColorsProxy const& o) = default;

    bool contains(size_t color_id) const {
        auto colorset = unitig.getData()->getUnitigColors(unitig);

        return colorset->contains(unitig, color_id);
    }

    UnitigColors::const_iterator begin() const {
        return unitig.getData()->getUnitigColors(unitig)->begin(unitig);
    }

    UnitigColors::const_iterator end() const {
        return unitig.getData()->getUnitigColors(unitig)->end();
    }

    size_t size() const {
        return unitig.getData()->getUnitigColors(unitig)->size(unitig);
    }

    size_t numKmersWithColor(size_t color_id) {
        return unitig.getData()->getUnitigColors(unitig)->size(unitig, color_id);
    }

private:
    PyfrostColoredUMap unitig;
};


void define_UnitigColors(py::module& m);

}

#endif //PYFROST_UNITIGCOLORS_H
