#ifndef PYFROST_ADJACENCYINNERDICT_H
#define PYFROST_ADJACENCYINNERDICT_H

#include <pybind11/pybind11.h>

#include "pyfrost.h"
#include "NodeIterator.h"

namespace py = pybind11;

namespace pyfrost {

enum class AdjacencyType {
    SUCCESSORS,
    PREDECESSORS
};

/**
 * This class emulates the inner dictionary of NetworkX's _adj dictionary, and contains
 * a specific node's successors or predecessors (depending on `type`).
 *
 * Accessing an element in this dictionary returns the edge data dictionary, which currently is hardcoded and
 * can't be modified. Bifrost's De Bruijn graph is node oriented and the edges don't have data associated.
 */
class AdjacencyInnerDict {
public:
    AdjacencyInnerDict(PyfrostCCDBG& _dbg, Kmer const& kmer, AdjacencyType _type);

private:
    PyfrostCCDBG& dbg;
    AdjacencyType type;
    PyfrostColoredUMap node;

public:
    using iterator_type = NodeIterator<decltype(node.getSuccessors().begin())>;

    PyfrostColoredUMap getUnitig(Kmer const& kmer);

    iterator_type begin() const;
    iterator_type end() const;

    bool contains(Kmer const& o) const;

    bool contains(const char* kmer) const {
        return contains(Kmer(kmer));
    }

    size_t numNeigbors() const;

    py::dict getEdgeDict(Kmer const& neighbor) const;

    py::dict getEdgeDict(char const* neighbor) const {
        return getEdgeDict(Kmer(neighbor));
    }
};

void define_AdjacencyInnerDict(py::module& m);

}

#endif //PYFROST_ADJACENCYINNERDICT_H
