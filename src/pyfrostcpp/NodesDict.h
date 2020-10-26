
#ifndef PYFROST_NODESDICT_H
#define PYFROST_NODESDICT_H

#include "pyfrost.h"
#include "NodeDataDict.h"
#include "NodeIterator.h"

namespace py = pybind11;

namespace pyfrost {

/**
 * This class represents NetworkX's _node attribute and can be used to access node metadata.
 */
class NodesDict {
public:
    explicit NodesDict(PyfrostCCDBG& dbg) : dbg(dbg) { }

    /**
     * This function searches for a unitig which starts with the given kmer.
     *
     * @param kmer
     */
    inline NodeDataDict findNode(Kmer const& kmer) {
        auto unitig = dbg.find(kmer, true).mappingToFullUnitig();

        if(unitig.isEmpty) {
            throw std::out_of_range("Node not found.");
        }

        return NodeDataDict(unitig);
    }

    inline NodeDataDict findNode(char const* kmer) {
        return findNode(Kmer(kmer));
    }

    /**
     * Return itself when calling with an existing kmer on unitig
     */
    bool contains(Kmer const& kmer) {
        auto unitig = dbg.find(kmer, true);

        return !unitig.isEmpty;
    }

    bool contains(char const* kmer) {
        return contains(Kmer(kmer));
    }

    /**
     * Iterate over unitigs (nodes)
     */
    inline auto begin() const {
        return NodeIterator<PyfrostCCDBG::iterator>(&dbg, dbg.begin(), true);
    }

    /**
     * End of the unitig (node) iterator
     */
    inline auto end() const {
        return NodeIterator<PyfrostCCDBG::iterator>(&dbg, dbg.end(), true);
    }

    /**
     * Iterate over unitigs (nodes) without reverse complements
     */
    inline auto begin_no_rc() const {
        return NodeIterator<PyfrostCCDBG::iterator>(&dbg, dbg.begin());
    }

    /**
     * End of the unitig (node) iterator without reverse complements
     */
    inline auto end_no_rc() const {
        return NodeIterator<PyfrostCCDBG::iterator>(&dbg, dbg.end());
    }

    size_t numNodes() const {
        // Also take into account reverse complements
        return dbg.size() * 2;
    }

private:
    PyfrostCCDBG& dbg;
};


void define_NodesDict(py::module& m);

}

#endif //PYFROST_NODESDICT_H
