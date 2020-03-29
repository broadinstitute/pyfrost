#include <pybind11/pybind11.h>
#include <CompactedDBG.hpp>

#include "Pyfrost.h"

#ifndef PYFROST_NODEVIEW_H
#define PYFROST_NODEVIEW_H

namespace py = pybind11;

namespace pyfrost {

class NodeView {
public:
    explicit NodeView(PyfrostCCDBG& dbg) : dbg(dbg) { }

    /**
     * This function searches for a unitig which starts with the given kmer.
     *
     * This function differs from `findKmer` because it doesn't consider internal k-mers.
     *
     * @param kmer
     */
    inline PyfrostColoredUMap findNode(Kmer const& kmer) {
        auto unitig = dbg.find(kmer, true);

        if(unitig.isEmpty) {
            throw std::out_of_range("Node not found.");
        }

        return unitig.mappingToFullUnitig();
    }

    inline PyfrostColoredUMap findNode(char const* kmer) {
        return findNode(Kmer(kmer));
    }

    /**
     * Return itself when calling with an existing kmer on unitig
     */
    inline PyfrostColoredUMap findNode(PyfrostColoredUMap const& unitig) {
        if(unitig.isEmpty || !(unitig.dist == 0 || unitig.dist == unitig.len)) {
            throw std::out_of_range("Node not found.");
        }

        return unitig.mappingToFullUnitig();
    }

    /**
     * Iterate over unitigs (nodes)
     */
    inline PyfrostCCDBG::iterator begin() const {
        return dbg.begin();
    }

    /**
     * End of the unitig (node) iterator
     */
    inline PyfrostCCDBG::iterator end() const {
        return dbg.end();
    }

    size_t numNodes() const {
        return dbg.size();
    }

private:
    PyfrostCCDBG& dbg;
};

void define_NodeView(py::module& m) {
    py::class_<NodeView>(m, "NodeView")
        // Access nodes with [] operator overloading
        .def("__getitem__", py::overload_cast<char const*>(&NodeView::findNode), py::is_operator())
        .def("__getitem__", py::overload_cast<Kmer const&>(&NodeView::findNode), py::is_operator())
        .def("__getitem__", py::overload_cast<PyfrostColoredUMap const&>(&NodeView::findNode), py::is_operator())

        .def("__len__", &NodeView::numNodes, py::is_operator())

        // Iterate over all nodes
        .def("__iter__", [](NodeView const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.begin(), self.end());
        }, py::keep_alive<0, 1>());
}


}

#endif //PYFROST_NODEVIEW_H
