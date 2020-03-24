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

    inline virtual decltype(unitig.getSuccessors().begin()) begin() const {
        return unitig.getSuccessors().begin();
    }

    inline virtual decltype(unitig.getSuccessors().end()) end() const {
        return unitig.getSuccessors().end();
    }
};

class PredecessorView : public AdjacencyViewBase {
public:
    explicit PredecessorView(PyfrostColoredUMap const& unitig) : AdjacencyViewBase(unitig) { }

    inline virtual decltype(unitig.getPredecessors().begin()) begin() const {
        return unitig.getPredecessors().begin();
    }

    inline virtual decltype(unitig.getPredecessors().end()) end() const {
        return unitig.getPredecessors().end();
    }
};


class AdjacencyProxy {
public:
    AdjacencyProxy(PyfrostCCDBG& dbg, AdjacencyType type) : dbg(dbg), type(type) { }

    inline AdjacencyViewBase* getView(PyfrostColoredUMap const& unitig) {
        if(type == AdjacencyType::SUCCESSORS) {
            return new SuccessorView(unitig);
        } else if(type == AdjacencyType::PREDECESSORS) {
            return new PredecessorView(unitig);
        } else {
            throw std::runtime_error("Unexpected value for AdjacencyProxy::type");
        }
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


void define_AdjacencyProxy(py::module& m) {
    py::enum_<AdjacencyType>(m, "AdjacencyType")
        .value("SUCCESSORS", AdjacencyType::SUCCESSORS)
        .value("PREDECESSORS", AdjacencyType::PREDECESSORS);

    py::class_<AdjacencyViewBase>(m, "AdjacencyViewBase")
        .def("__len__", &AdjacencyViewBase::numNeigbors, py::is_operator())

        .def("__contains__", py::overload_cast<char const*>(&AdjacencyViewBase::contains), py::is_operator())
        .def("__contains__", py::overload_cast<Kmer const&>(&AdjacencyViewBase::contains), py::is_operator())
        .def("__contains__", py::overload_cast<PyfrostColoredUMap const&>(&AdjacencyViewBase::contains), py::is_operator())

        .def("__getitem__", &AdjacencyViewBase::getEdgeDict);

    py::class_<SuccessorView, AdjacencyViewBase>(m, "SuccessorView")
        .def(py::init<PyfrostColoredUMap const&>())
        .def("__iter__", [] (SuccessorView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    py::class_<PredecessorView, AdjacencyViewBase>(m, "PredecessorView")
        .def(py::init<PyfrostColoredUMap const&>())
        .def("__iter__", [] (PredecessorView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    py::class_<AdjacencyProxy>(m, "AdjacencyProxy")
        .def("__len__", &AdjacencyProxy::numNodes, py::is_operator())

        .def("__contains__", py::overload_cast<PyfrostColoredUMap const&>(&AdjacencyProxy::contains))
        .def("__contains__", py::overload_cast<Kmer const&>(&AdjacencyProxy::contains))
        .def("__contains__", py::overload_cast<char const*>(&AdjacencyProxy::contains))

        .def("__getitem__", &AdjacencyProxy::getView, py::is_operator())

        .def("__iter__", [] (AdjacencyProxy const& self) {
            return py::make_iterator(self.begin(), self.end(), py::return_value_policy::copy);
        }, py::keep_alive<0, 1>());
}

}



#endif //PYFROST_ADJACENCYVIEW_H
