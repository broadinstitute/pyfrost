#include <pybind11/pybind11.h>
#include <NeighborIterator.hpp>

#include "Pyfrost.h"

#ifndef PYFROST_ADJACENCYVIEW_H
#define PYFROST_ADJACENCYVIEW_H

namespace py = pybind11;

namespace pyfrost {

enum class AdjacencyType {
    SUCCESSORS,
    PREDECESSORS
};

class AdjacencyViewBase {
public:
    explicit AdjacencyViewBase(PyfrostColoredUMap const& unitig) : unitig(unitig) {}
    virtual ~AdjacencyViewBase() { }

protected:
    PyfrostColoredUMap unitig;

public:
    inline virtual decltype(unitig.getSuccessors().begin()) begin() const = 0;
    inline virtual decltype(unitig.getSuccessors().end()) end() const = 0;

    bool contains(PyfrostColoredUMap const& o) {
        // We have at most 4 neighbors so this should be quick enough, and could be considered constant time.
        for(auto const& neighbor : *this) {
            if(o == neighbor) {
                return true;
            }
        }

        return false;
    }

    bool contains(Kmer const& kmer) {
        auto o = unitig.getGraph()->find(kmer, true);
        return contains(o);
    }

    bool contains(const char* kmer) {
        auto o = unitig.getGraph()->find(Kmer(kmer), true);
        return contains(o);
    }

    py::dict getEdgeDict(PyfrostColoredUMap const& o) {
        return py::dict();
    }

    size_t numNeigbors() const {
        size_t num = 0;
        for(auto const& n : *this) {
            ++num;
        }

        return num;
    }

};

class SuccessorView : public AdjacencyViewBase {
public:
    explicit SuccessorView(PyfrostColoredUMap const& unitig) : AdjacencyViewBase(unitig) { }

    inline decltype(unitig.getSuccessors().begin()) begin() const override {
        return unitig.getSuccessors().begin();
    }

    inline decltype(unitig.getSuccessors().end()) end() const override {
        return unitig.getSuccessors().end();
    }
};

class PredecessorView : public AdjacencyViewBase {
public:
    explicit PredecessorView(PyfrostColoredUMap const& unitig) : AdjacencyViewBase(unitig) { }

    inline decltype(unitig.getPredecessors().begin()) begin() const override {
        return unitig.getPredecessors().begin();
    }

    inline decltype(unitig.getPredecessors().end()) end() const override {
        return unitig.getPredecessors().end();
    }
};


class AdjacencyProxy {
public:
    AdjacencyProxy(PyfrostCCDBG& dbg, AdjacencyType type) : dbg(dbg), type(type) { }

    inline AdjacencyViewBase* getView(PyfrostColoredUMap const& unitig) {
        if(unitig.isEmpty) {
            throw std::out_of_range("Node does not exist in the graph.");
        }

        if(type == AdjacencyType::SUCCESSORS) {
            return new SuccessorView(unitig);
        } else {
            return new PredecessorView(unitig);
        }
    }

    inline AdjacencyViewBase* getView(Kmer const& kmer) {
        auto unitig = dbg.find(kmer);
        return getView(unitig);
    }

    inline AdjacencyViewBase* getView(char const* kmer) {
        auto unitig = dbg.find(Kmer(kmer));
        return getView(unitig);
    }

    inline auto begin() const {
        return dbg.begin();
    }

    inline auto end() const {
        return dbg.end();
    }

    inline bool contains(PyfrostColoredUMap const& unitig) {
        return unitig.dist == 0 && !unitig.isEmpty;
    }

    inline bool contains(Kmer const& kmer) {
        auto um = dbg.find(kmer, true);
        return contains(um);
    }

    inline bool contains(char const* kmer) {
        auto um = dbg.find(Kmer(kmer), true);
        return contains(um);
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
