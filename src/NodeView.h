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

    inline PyfrostColoredUMap find_kmer(Kmer const& kmer) {
        return dbg.find(kmer, false);
    }

    inline PyfrostColoredUMap find_kmer(char const* kmer) {
        Kmer km = Kmer(kmer);
        return dbg.find(km, false);
    }

    /**
     * This function searches for a unitig which starts with the given kmer.
     *
     * This function differs from `find_kmer` because it doesn't consider internal k-mers.
     *
     * @param kmer
     */
    inline PyfrostColoredUMap find_node(Kmer const& kmer) {
        return dbg.find(kmer, true);
    }

    inline PyfrostColoredUMap find_node(char const* kmer) {
        Kmer km = Kmer(kmer);
        return dbg.find(km, true);
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

private:
    PyfrostCCDBG& dbg;
};

void define_NodeView(py::module& m) {
    py::class_<NodeView>(m, "NodeView")
        // Access nodes with [] operator overloading
        .def("__getitem__", py::overload_cast<char const*>(&NodeView::find_node), py::is_operator())
        .def("__getitem__", py::overload_cast<Kmer const&>(&NodeView::find_node), py::is_operator())

        // Iterate over all nodes
        .def("__iter__", [](NodeView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());
}


}

#endif //PYFROST_NODEVIEW_H
