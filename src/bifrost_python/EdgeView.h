#ifndef PYFROST_EDGEVIEW_H
#define PYFROST_EDGEVIEW_H

#include <numeric>

#include "Pyfrost.h"
#include "NodeBunchIter.h"

namespace pyfrost {

/**
 * Iterate over all outgoing edges in the graph, optionally with edge metadata
 *
 * \tparam T Node iterator type, can be used to iterate over different node containers. Dereference of this iterator
 *           type should return a UnitigMap
 */
template<typename T>
class OutEdgeIterator {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = py::object;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    explicit OutEdgeIterator(T const& node) : node_iter(node),
        edge_iter(node->getSuccessors().begin()),
        edge_end(node->getSuccessors().end()) { }
    OutEdgeIterator(OutEdgeIterator const& o) = default;
    OutEdgeIterator(OutEdgeIterator&& o) = default;

    OutEdgeIterator& operator++();
    OutEdgeIterator operator++(int);

    value_type operator*();

    bool operator==(OutEdgeIterator const& o) const;
    bool operator!=(OutEdgeIterator const& o) const;

private:
    T node_iter;
    PyfrostColoredUMap::neighbor_iterator edge_iter;
    PyfrostColoredUMap::neighbor_iterator edge_end;
};

template<typename T>
OutEdgeIterator<T>& OutEdgeIterator<T>::operator++() {
    do {
        ++edge_iter;

        // All edges of current node visited, move to next node
        if(edge_iter == edge_end) {
            ++node_iter;

            edge_iter = node_iter->getSuccessors().begin();
            edge_end = node_iter->getSuccessors().end();
        }
    } while (!node_iter->isEmpty && edge_iter == edge_end);

    return *this;
}

template<typename T>
OutEdgeIterator<T> OutEdgeIterator<T>::operator++(int) {
    OutEdgeIterator tmp(*this);

    operator++();
    return tmp;
}

template<typename T>
typename OutEdgeIterator<T>::value_type OutEdgeIterator<T>::operator*() {
    return py::make_tuple(*node_iter, *edge_iter);
}

template<typename T>
bool OutEdgeIterator<T>::operator==(OutEdgeIterator<T> const& o) const {
    return node_iter == o.node_iter && edge_iter == o.edge_iter;
}

template<typename T>
bool OutEdgeIterator<T>::operator!=(OutEdgeIterator<T> const& o) const {
    return !operator==(o);
}


class BifrostDiGraph;

/**
 * Simulate NetworkX's edge view
 *
 * Pyfrost, however, doesn't support user metadata on edges, and will return a hardcoded dict with some basic metadata.
 */
class EdgeView {
public:
    explicit EdgeView(BifrostDiGraph& dbg);
    EdgeView(EdgeView const& o) = default;
    EdgeView(EdgeView&& o) = default;

    py::dict get(py::tuple const& edge);
    py::dict get(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2);

    OutEdgeIterator<PyfrostCCDBG::iterator> begin() const;
    OutEdgeIterator<PyfrostCCDBG::iterator> end() const;

    size_t size() const;

private:
    BifrostDiGraph& dbg;
};


void define_EdgeView(py::module& m);

}

#endif //PYFROST_EDGEVIEW_H
