#ifndef PYFROST_ADJACENCYVIEW_H
#define PYFROST_ADJACENCYVIEW_H

#include <pybind11/pybind11.h>
#include <NeighborIterator.hpp>

#include "Pyfrost.h"
#include "NodeIterator.h"
#include "EdgeView.h"

namespace py = pybind11;

namespace pyfrost {

enum class AdjacencyType {
    SUCCESSORS,
    PREDECESSORS
};

class AdjacencyViewBase {
public:
    AdjacencyViewBase(PyfrostCCDBG& dbg, Kmer const& kmer) : dbg(dbg) {
        node = getUnitig(kmer);
    }
    virtual ~AdjacencyViewBase() { }

protected:
    PyfrostCCDBG& dbg;
    PyfrostColoredUMap node;

public:
    using iterator_type = NodeIterator<decltype(node.getSuccessors().begin())>;

    virtual iterator_type begin() const = 0;
    virtual iterator_type end() const = 0;

    PyfrostColoredUMap getUnitig(Kmer const& kmer) {
        auto unitig = dbg.find(kmer, true).mappingToFullUnitig();

        if(unitig.isEmpty) {
            throw std::out_of_range("Node does not exist in the graph.");
        }

        return unitig;
    }

    bool contains(Kmer const& o) {
        // We have at most 4 neighbors so this should be quick enough, and could be considered constant time.
        for(Kmer const& neighbor : *this) {
            if(o == neighbor && !is_kmer_empty(neighbor)) {
                return true;
            }
        }

        return false;
    }

    bool contains(const char* kmer) {
        return contains(Kmer(kmer));
    }

    size_t numNeigbors() const {
        size_t num = 0;
        for(auto const& n : *this) {
            ++num;
        }

        return num;
    }

    virtual py::dict getEdgeDict(Kmer const& neighbor) const = 0;

    py::dict getEdgeDict(char const* neighbor) const {
        return getEdgeDict(Kmer(neighbor));
    }
};

class SuccessorView : public AdjacencyViewBase {
public:
    explicit SuccessorView(PyfrostCCDBG& dbg, Kmer const& kmer) : AdjacencyViewBase(dbg, kmer) { }

    iterator_type begin() const override {
        return {&dbg, node.getSuccessors().begin(), false};
    }

    iterator_type end() const override {
        return {&dbg, node.getSuccessors().end(), false};
    }

    py::dict getEdgeDict(Kmer const& neighbor) const override {
        return makeEdgeDataDict(dbg, node.getMappedHead(), neighbor);
    }
};

class PredecessorView : public AdjacencyViewBase {
public:
    explicit PredecessorView(PyfrostCCDBG& dbg, Kmer const& kmer) : AdjacencyViewBase(dbg, kmer) { }

    iterator_type begin() const override {
        return {&dbg, node.getPredecessors().begin(), false};
    }

    iterator_type end() const override {
        return {&dbg, node.getPredecessors().end(), false};
    }

    py::dict getEdgeDict(Kmer const& neighbor) const override {
        return makeEdgeDataDict(dbg, neighbor, node.getMappedHead());
    }
};


class AdjacencyProxy {
public:
    AdjacencyProxy(PyfrostCCDBG& dbg, AdjacencyType type) : dbg(dbg), type(type) { }

    inline AdjacencyViewBase* getView(Kmer const& kmer) {
        if(type == AdjacencyType::SUCCESSORS) {
            return new SuccessorView(dbg, kmer);
        } else {
            return new PredecessorView(dbg, kmer);
        }
    }

    inline AdjacencyViewBase* getView(char const* kmer) {
        return getView(Kmer(kmer));
    }

    auto begin() const {
        return NodeIterator<PyfrostCCDBG::iterator>(&dbg, dbg.begin(), false);
    }

    auto end() const {
        return NodeIterator<PyfrostCCDBG::iterator>(&dbg, dbg.end(), false);
    }

    inline bool contains(Kmer const& kmer) {
        auto um = dbg.find(kmer, true);
        return !um.isEmpty;
    }

    inline bool contains(char const* kmer) {
        return contains(Kmer(kmer));
    }

    size_t numNodes() const {
        return dbg.size();
    }

private:
    PyfrostCCDBG& dbg;
    AdjacencyType type;
};


void define_AdjacencyProxy(py::module& m);

}

#endif //PYFROST_ADJACENCYVIEW_H
