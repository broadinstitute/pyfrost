#ifndef PYFROST_UNITIGCOVERAGE_H
#define PYFROST_UNITIGCOVERAGE_H

#include <limits>

#include "Pyfrost.h"

namespace pyfrost {

class UnitigCoverageIterator {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = size_t;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    explicit UnitigCoverageIterator(PyfrostColoredUMap& unitig)
        : unitig(unitig), pos(std::numeric_limits<size_t>::max())
    {
    }

    UnitigCoverageIterator(PyfrostColoredUMap& unitig, size_t pos) : unitig(unitig), pos(pos)
    {
    }

    UnitigCoverageIterator(UnitigCoverageIterator const& o) = default;

    UnitigCoverageIterator(UnitigCoverageIterator&& o) = default;

    value_type operator*()
    {
        if(unitig.isEmpty) {
            return 0;
        }

        return unitig.getCoverage(pos);
    }

    void operator++()
    {
        ++pos;
    }

    UnitigCoverageIterator operator++(int)
    {
        UnitigCoverageIterator tmp(*this);

        operator++();
        return tmp;
    }

    bool operator==(UnitigCoverageIterator const& o)
    {
        return pos == o.pos;
    }

    bool operator!=(UnitigCoverageIterator const& o)
    {
        return pos != o.pos;
    }

private:
    PyfrostColoredUMap& unitig;
    size_t pos;
};

class UnitigCoverageProxy {
public:
    explicit UnitigCoverageProxy(PyfrostColoredUMap const& src) : unitig(src)
    {

    }

    size_t size() const
    {
        return unitig.size - unitig.getGraph()->getK() + 1;
    }

    size_t getCoverage(long pos)
    {
        if(unitig.isEmpty) {
            return 0;
        }

        size_t new_pos = pos < 0 ? unitig.size - unitig.getGraph()->getK() + 1 + pos : static_cast<size_t>(pos);
        return unitig.getCoverage(new_pos);
    }

    UnitigCoverageIterator begin()
    {
        return UnitigCoverageIterator(unitig, 0);
    }

    UnitigCoverageIterator end()
    {
        return UnitigCoverageIterator(unitig, size());
    }

private:
    PyfrostColoredUMap unitig;
};


void define_UnitigCoverageProxy(py::module& m);

}


#endif //PYFROST_UNITIGCOVERAGE_H
