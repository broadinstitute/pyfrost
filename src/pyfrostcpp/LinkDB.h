#ifndef PYFROST_LINKDB_H
#define PYFROST_LINKDB_H

#include <string>
#include <thread>
#include <queue>
#include <robin_hood.h>
#include <tl/optional.hpp>

#include "pyfrost.h"
#include "JunctionTree.h"

namespace pyfrost {

using junction_tree_map = robin_hood::unordered_map<Kmer, unique_ptr<JunctionTreeNode>>;

class LinkDB {
public:
    LinkDB() = default;
    explicit LinkDB(size_t _color) : color(_color) { }

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
    virtual JunctionTreeNode& getLinks(Kmer const& kmer) = 0;

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
    virtual JunctionTreeNode& createOrGetTree(Kmer const& kmer) = 0;

    /**
     * Merge links from another database with this link database
     */
    virtual void merge(LinkDB& other) {
        for(auto& pair : other.getJunctionTrees()) {
            if(!pair.second->getChildren().empty()) {
                auto& tree = createOrGetTree(pair.first);
                tree.merge(*pair.second);
            }
        }
    }

    /**
     * Get the color associated with this link database.
     */
    tl::optional<size_t> getColor() const {
        return color;
    }

    /**
     * Set the color associated with this link database.
     */
    void setColor(size_t _color) {
        color = _color;
    }

    /**
     * To be called after loading from a file, which doesn't store the parent pointers.
     */
    void fixTreeParents() {
        for(auto& tree : getJunctionTrees()) {
            tree.second->fixParents();
        }
    }

protected:
    /// With which color is this link database associated?
    tl::optional<size_t> color;

};

void define_LinkDB(py::module& m);

}

#endif //PYFROST_LINKDB_H
