#ifndef PYFROST_LINKDB_H
#define PYFROST_LINKDB_H

#include <string>
#include <thread>
#include <queue>
#include <robin_hood.h>

#include "pyfrost.h"
#include "JunctionTree.h"

namespace pyfrost {

using junction_tree_map = robin_hood::unordered_map<Kmer, std::shared_ptr<JunctionTreeNode>>;

class LinkDB {
public:
    LinkDB() = default;

    LinkDB(LinkDB const& o) = default;
    LinkDB(LinkDB&& o) noexcept = default;

    virtual ~LinkDB() = default;

    LinkDB& operator=(LinkDB const& o) = default;

    /**
     * Check if the database has links for a given kmer
     */
    virtual bool hasLinks(Kmer const& kmer) const = 0;

    /**
     * Get a JunctionTreeNode representing all links for a given k-mer
     */
    virtual std::shared_ptr<JunctionTreeNode> getLinks(Kmer const& kmer) = 0;

    /**
     * Total number of JunctionTrees in this database
     */
    virtual size_t numTrees() const = 0;

    /**
     * Return all <kmer, JunctionTreeNode> instances in this database.
     */
    virtual junction_tree_map& getJunctionTrees() = 0;

    /**
     * Create a new junction tree for the given k-mer in the database, or return the existing one if available.
     */
    virtual std::shared_ptr<JunctionTreeNode> createOrGetTree(Kmer const& kmer) = 0;

};

void define_LinkDB(py::module& m);

}

#endif //PYFROST_LINKDB_H
