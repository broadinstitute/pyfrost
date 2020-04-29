#ifndef PYFROST_NODEBUNCHITER_H
#define PYFROST_NODEBUNCHITER_H

#include "Pyfrost.h"

namespace pyfrost {

class BifrostDiGraph;

class NodeBunchIter {
private:
    py::iterator wrapped;

public:
    using iterator_category = std::input_iterator_tag;
    using value_type = PyfrostColoredUMap;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type const*;

    NodeBunchIter();
    explicit NodeBunchIter(py::iterator iter);
    NodeBunchIter(NodeBunchIter const& o) = default;
    NodeBunchIter(NodeBunchIter&& o) = default;

    value_type operator*();
    pointer operator->() const;

    NodeBunchIter& operator++();
    NodeBunchIter operator++(int);

    bool operator==(NodeBunchIter const& o) const;
    bool operator!=(NodeBunchIter const& o) const;
};

}

#endif //PYFROST_NODEBUNCHITER_H
