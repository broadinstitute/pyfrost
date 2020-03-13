#include <CompactedDBG.hpp>

#include "Pyfrost.h"
#include "NodeView.h"

#ifndef PYFROST_BIFROSTGRAPH_H
#define PYFROST_BIFROSTGRAPH_H

namespace py = pybind11;

namespace pyfrost {

class BifrostDiGraph {
public:
    BifrostDiGraph() : nodes(dbg);
    BifrostDiGraph(BifrostDiGraph const& o) : dbg(o.dbg), nodes(dbg) { }
    BifrostDiGraph(BifrostDiGraph&& o) : dbg(std::move(o.dbg)), nodes(dbg) { }
    BifrostDiGraph(PyfrostCCDBG&& dbg) : dbg(std::move(dbg)), nodes(dbg) { }

    NodeView& getNodeView();

private:
    PyfrostCCDBG dbg;
    NodeView nodes;

};

void define_BifrostDiGraph(py::module& m) {
    py::class_<BifrostDiGraph>(m, "BifrostDiGraph")
        .def(py::init<>())
        .def(py::init<BifrostDiGraph const&>());

}

}


#endif //PYFROST_BIFROSTGRAPH_H
