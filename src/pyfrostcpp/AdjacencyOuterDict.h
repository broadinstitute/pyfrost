#ifndef PYFROST_ADJACENCYVIEW_H
#define PYFROST_ADJACENCYVIEW_H

#include <pybind11/pybind11.h>
#include <NeighborIterator.hpp>

#include "pyfrost.h"
#include "AdjacencyInnerDict.h"
#include "NodeIterator.h"

namespace py = pybind11;

namespace pyfrost {


/**
 * This class emulates the outer dict of NetworkX's _adj attribute, and represents a dictionary with all
 * successors/predecessors keyed by node.
 *
 * Accessing an element in this dictionary returns an AdjacencyInnerDict, which is another dictionary that lists
 * successors/predecessors for a specific node.
 */
class AdjacencyOuterDict {
public:
    AdjacencyOuterDict(PyfrostCCDBG& _dbg, AdjacencyType _type) : dbg(_dbg), type(_type) { };

    AdjacencyInnerDict getInnerDict(Kmer const& kmer) {
        return {dbg, kmer, type};
    }

    AdjacencyInnerDict getInnerDict(char const* kmer) {
        return getInnerDict(Kmer(kmer));
    }

    auto begin() const {
        return NodeIterator<PyfrostCCDBG::iterator>(&dbg, dbg.begin(), false);
    }

    auto end() const {
        return NodeIterator<PyfrostCCDBG::iterator>(&dbg, dbg.end(), false);
    }

    bool contains(Kmer const& kmer) {
        auto um = dbg.find(kmer, true);
        return !um.isEmpty;
    }

    bool contains(char const* kmer) {
        return contains(Kmer(kmer));
    }

    size_t numNodes() const {
        return dbg.size();
    }

private:
    PyfrostCCDBG& dbg;
    AdjacencyType type;
};


void define_AdjacencyOuterDict(py::module& m);

}

#endif //PYFROST_ADJACENCYVIEW_H
