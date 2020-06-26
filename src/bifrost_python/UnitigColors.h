#ifndef PYFROST_UNITIGCOLORS_H
#define PYFROST_UNITIGCOLORS_H

#include "Pyfrost.h"

namespace py = pybind11;

namespace pyfrost {

class UnitigColorIterator {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = size_t;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type const*;

    UnitigColorIterator();
    UnitigColorIterator(UnitigColors::const_iterator const& unitig);
    UnitigColorIterator(UnitigColorIterator const& o) = default;
    UnitigColorIterator(UnitigColorIterator&&) = default;

    value_type operator*() const;
    pointer operator->() const;

    void operator++();
    UnitigColorIterator operator++(int);

    bool operator==(UnitigColorIterator const& o);
    bool operator!=(UnitigColorIterator const& o);

private:
    UnitigColors::const_iterator iter;
    value_type current;
};

class UnitigColorsProxy {
public:
    explicit UnitigColorsProxy(PyfrostColoredUMap const& unitig) : unitig(unitig) { }

    UnitigColorsProxy(UnitigColorsProxy const& o) = default;
    UnitigColorsProxy(UnitigColorsProxy&& o) = default;

    UnitigColorsProxy& operator=(UnitigColorsProxy const& o) = default;

    UnitigColorsProxy getColorsAtPos(int pos) {
        size_t new_pos = pos < 0 ? unitig.size - unitig.getGraph()->getK() + 1 + pos : static_cast<size_t>(pos);
        auto new_unitig = unitig.getKmerMapping(new_pos);

        if(new_unitig.isEmpty) {
            throw py::key_error("Index out of range");
        }

        return UnitigColorsProxy(new_unitig);
    }

    bool contains(size_t color_id) const {
        auto colorset = unitig.getData()->getUnitigColors(unitig);

        return colorset->contains(unitig, color_id);
    }

    UnitigColorIterator begin() const {
        return UnitigColorIterator(unitig.getData()->getUnitigColors(unitig)->begin(unitig));
    }

    UnitigColorIterator end() const {
        return UnitigColorIterator(unitig.getData()->getUnitigColors(unitig)->end());
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
