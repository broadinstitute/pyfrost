#ifndef PYFROST_MEMLINKDB_H
#define PYFROST_MEMLINKDB_H

#include "pyfrost.h"
#include "JunctionTree.h"
#include "LinkDB.h"

namespace pyfrost {

class MemLinkDB : public LinkDB {
public:
    MemLinkDB() { }

    MemLinkDB(MemLinkDB const& o) = default;
    MemLinkDB(MemLinkDB&& o) noexcept = default;

    MemLinkDB& operator=(MemLinkDB const& o) = default;

    std::shared_ptr<JunctionTreeNode> getLinks(Kmer const& kmer) override;

    size_t numTrees() const override {
        return junction_trees.size();
    }

    junction_tree_map& getJunctionTrees() override {
        return junction_trees;
    }

    std::shared_ptr<JunctionTreeNode> createOrGetTree(Kmer const& kmer) override;

private:
    /// A map which stores the junction trees associated with each k-mer
    junction_tree_map junction_trees;
};

void define_MemLinkDB(py::module& m);

}

#endif //PYFROST_MEMLINKDB_H
