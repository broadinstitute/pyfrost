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

class AdjacencyView {

public:
    explicit AdjacencyView(PyfrostColoredUMap const& unitig, AdjacencyType type) : unitig(unitig) {}

    inline auto begin() const {
        if(type == AdjacencyType::SUCCESSORS) {
            return unitig.getSuccessors().begin();
        }
        else {
            return unitig.getPredecessors().begin();
        }
    }

    inline auto end() const {
        if(type == AdjacencyType::SUCCESSORS) {
            return unitig.getSuccessors().end();
        } else {
            return unitig.getPredecessors().end();
        }
    }

private:
    PyfrostColoredUMap const& unitig;
    AdjacencyType type;

};


class AdjacencyProxy {
public:
    explicit AdjacencyProxy(PyfrostCCDBG& dbg) : dbg(dbg), type(AdjacencyType::SUCCESSORS) { }
    AdjacencyProxy(PyfrostCCDBG& dbg, AdjacencyType type) : dbg(dbg), type(type) { }

    inline AdjacencyView* getView(PyfrostColoredUMap const& unitig) {
        return new AdjacencyView(unitig, type);
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
        return this->contains(um);
    }

    inline bool contains(char const* kmer) {
        Kmer km(kmer);
        auto um = dbg.find(km, true);
        return this->contains(um);
    }

private:
    PyfrostCCDBG& dbg;
    AdjacencyType type;
};

void define_AdjacencyProxy(py::module& m) {
    py::class_<AdjacencyView>(m, "AdjacencyView")
        .def("__iter__", [] (AdjacencyView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    py::class_<AdjacencyProxy>(m, "AdjacencyProxy")
        .def("__getitem__", &AdjacencyProxy::getView, py::is_operator())
        .def("__iter__", [] (AdjacencyProxy const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        .def("__contains__", py::overload_cast<PyfrostColoredUMap const&>(&AdjacencyProxy::contains))
        .def("__contains__", py::overload_cast<Kmer const&>(&AdjacencyProxy::contains))
        .def("__contains__", py::overload_cast<char const*>(&AdjacencyProxy::contains));
}

}



#endif //PYFROST_ADJACENCYVIEW_H
