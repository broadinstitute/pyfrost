#ifndef PYFROST_NODEBUNCHITER_H
#define PYFROST_NODEBUNCHITER_H

#include "Pyfrost.h"

namespace pyfrost {

class BifrostDiGraph;

class NodeBunchIter {
private:
    using wrapped_iter_t = py::iterator;

    BifrostDiGraph const &graph;
    wrapped_iter_t wrapped;

public:
    using iterator_category = std::input_iterator_tag;
    using value_type = PyfrostColoredUMap;
    using difference_type = std::ptrdiff_t;
    using reference = value_type &;
    using pointer = value_type *;

    NodeBunchIter(BifrostDiGraph const &graph, py::iterator iter) : graph(graph), wrapped(std::move(iter)) {}
    NodeBunchIter(NodeBunchIter const &o) = default;
    NodeBunchIter(NodeBunchIter &&o) = default;

    value_type operator*();

    NodeBunchIter &operator++();
    NodeBunchIter operator++(int);

    bool operator==(NodeBunchIter const &o);
    bool operator!=(NodeBunchIter const &o);
};

}

#endif //PYFROST_NODEBUNCHITER_H
