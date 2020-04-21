#ifndef PYFROST_EDGEVIEW_H
#define PYFROST_EDGEVIEW_H

#include <numeric>

#include "Pyfrost.h"
#include "NodeBunchIter.h"

namespace pyfrost {

py::dict makeEdgeDataDict(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2);

/**
 * Iterate over all outgoing out_edges in the graph, optionally with edge metadata
 *
 * \tparam T Node iterator type, can be used to iterate over different node containers. Dereference of this iterator
 *           type should return a UnitigMap
 */
template<typename T, bool outgoing=true>
class EdgeIterator {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = py::object;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    EdgeIterator() = default;

    EdgeIterator(T const& begin, T const& end) : node_iter(begin), node_end(end) {
        // find a node with outgoing or incoming edges
        while(node_iter != node_end) {
            size_t num = (outgoing ? node_iter->getSuccessors().cardinality() :
                          node_iter->getPredecessors().cardinality());

            if(num > 0) {
                edge_iter = outgoing ? node_iter->getSuccessors().begin() : node_iter->getPredecessors().begin();
                edge_end = outgoing ? node_iter->getSuccessors().end() : node_iter->getPredecessors().end();

                break;
            } else {
                ++node_iter;
            }
        }
    }

    EdgeIterator(EdgeIterator const& o) = default;
    EdgeIterator(EdgeIterator&& o) = default;

    EdgeIterator& operator++();
    EdgeIterator operator++(int);

    value_type operator*();

    bool operator==(EdgeIterator const& o) const;
    bool operator!=(EdgeIterator const& o) const;

private:
    T node_iter;
    T node_end;
    PyfrostColoredUMap::neighbor_iterator edge_iter;
    PyfrostColoredUMap::neighbor_iterator edge_end;
};

template<typename T, bool outgoing>
EdgeIterator<T, outgoing>& EdgeIterator<T, outgoing>::operator++() {
    if(node_iter->isEmpty || node_iter == node_end) {
        return *this;
    }

    do {
        // Move to next edge if possible
        ++edge_iter;

        if(edge_iter == edge_end) {
            // All edges of current node visited, move to next node
            ++node_iter;

            if(node_iter == node_end) {
                break;
            }

            if(outgoing) {
                edge_iter = node_iter->getSuccessors().begin();
                edge_end = node_iter->getSuccessors().end();
            } else {
                edge_iter = node_iter->getPredecessors().begin();
                edge_end = node_iter->getPredecessors().end();
            }
        }
    } while (edge_iter == edge_end);

    return *this;
}

template<typename T, bool outgoing>
EdgeIterator<T, outgoing> EdgeIterator<T, outgoing>::operator++(int) {
    EdgeIterator<T, outgoing> tmp(*this);

    operator++();
    return tmp;
}

template<typename T, bool outgoing>
typename EdgeIterator<T, outgoing>::value_type EdgeIterator<T, outgoing>::operator*() {
    return py::make_tuple(*node_iter, *edge_iter);
}

template<typename T, bool outgoing>
bool EdgeIterator<T, outgoing>::operator==(EdgeIterator<T, outgoing> const& o) const {
    if(node_iter == node_end && o.node_iter == o.node_end) {
        return true;
    } else if(node_iter->isEmpty && o.node_iter->isEmpty) {
        return true;
    } else {
        return node_iter == o.node_iter && edge_iter == o.edge_iter;
    }
}

template<typename T, bool outgoing>
bool EdgeIterator<T, outgoing>::operator!=(EdgeIterator<T, outgoing> const& o) const {
    return !operator==(o);
}


class BifrostDiGraph;

class OutEdgeView {
public:
    explicit OutEdgeView(BifrostDiGraph& dbg);
    OutEdgeView(OutEdgeView const& o) = default;
    OutEdgeView(OutEdgeView&& o) = default;

    py::dict get(py::tuple const& edge);
    py::dict get(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2);

    EdgeIterator<PyfrostCCDBG::iterator, true> begin() const;
    EdgeIterator<PyfrostCCDBG::iterator, true> end() const;

    bool contains(py::tuple const& edge);

    virtual size_t size() const;

private:
    BifrostDiGraph& dbg;
};


class InEdgeView {
public:
    explicit InEdgeView(BifrostDiGraph& dbg);
    InEdgeView(InEdgeView const& o) = default;
    InEdgeView(InEdgeView&& o) = default;

    py::dict get(py::tuple const& edge);
    py::dict get(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2);

    EdgeIterator<PyfrostCCDBG::iterator, false> begin() const;
    EdgeIterator<PyfrostCCDBG::iterator, false> end() const;

    bool contains(py::tuple const& edge);

    size_t size() const;

private:
    BifrostDiGraph& dbg;
};

void define_EdgeView(py::module& m);

}

#endif //PYFROST_EDGEVIEW_H
