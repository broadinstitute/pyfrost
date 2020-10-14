#ifndef PYFROST_UNITIGCOLORS_H
#define PYFROST_UNITIGCOLORS_H

#include "pyfrost.h"

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

    UnitigColorIterator& operator++();
    UnitigColorIterator operator++(int);

    bool operator==(UnitigColorIterator const& o);
    bool operator!=(UnitigColorIterator const& o);

private:
    UnitigColors::const_iterator iter;
    value_type current;
};

class UnitigColorsProxy {
public:
    explicit UnitigColorsProxy(PyfrostColoredUMap const& _unitig, bool _rev_compl = false) :
        unitig(_unitig), rev_compl(_rev_compl), colorset(unitig.getData()->getUnitigColors(unitig))
    {
        if(colorset == nullptr) {
            throw std::runtime_error("Invalid colorset for unitig! Got a nullptr...");
        }
    }

    UnitigColorsProxy(UnitigColorsProxy const& o) = default;
    UnitigColorsProxy(UnitigColorsProxy&& o) = default;

    UnitigColorsProxy getColorsAtPos(long pos) {
        size_t new_pos;

        if(rev_compl) {
            new_pos = pos >= 0 ? unitig.size - unitig.getGraph()->getK() + 1 + pos : static_cast<size_t>(pos);
        } else {
            new_pos = pos < 0 ? unitig.size - unitig.getGraph()->getK() + 1 + pos : static_cast<size_t>(pos);
        }

        auto new_unitig = unitig.getKmerMapping(new_pos);

        if(new_unitig.isEmpty) {
            throw py::key_error("Index out of range");
        }

        return UnitigColorsProxy(new_unitig);
    }

    bool contains(size_t color_id) const {
        return colorset->contains(unitig, color_id);
    }

    void add(size_t color_id) {
        colorset->add(unitig, color_id);
    }

    void discard(size_t color_id) {
        colorset->remove(unitig, color_id);
    }

    UnitigColorIterator begin() const {
        return UnitigColorIterator(colorset->begin(unitig));
    }

    UnitigColorIterator end() const {
        return UnitigColorIterator(colorset->end());
    }

    size_t size() const {
        return colorset->size(unitig);
    }

    size_t numKmersWithColor(size_t color_id) {
        return colorset->size(unitig, color_id);
    }

private:
    PyfrostColoredUMap unitig;
    bool rev_compl;
    UnitigColors* colorset;
};


void define_UnitigColors(py::module& m);

}

#endif //PYFROST_UNITIGCOLORS_H
