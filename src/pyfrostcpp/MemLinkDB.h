#ifndef PYFROST_MEMLINKDB_H
#define PYFROST_MEMLINKDB_H

#include "pyfrost.h"
#include "JunctionTree.h"
#include "LinkDB.h"
#include "Serialize.h"

namespace pyfrost {

class MemLinkDB : public LinkDB {
public:
    MemLinkDB() : LinkDB() { }
    explicit MemLinkDB(size_t _color) : LinkDB(_color) { }
    MemLinkDB(MemLinkDB&& o) noexcept = default;

    MemLinkDB& operator=(MemLinkDB const& o) = default;

    bool hasLinks(Kmer const& kmer) const override {
        return junction_trees.contains(kmer);
    }

    JunctionTreeNode& getLinks(Kmer const& kmer) override;

    size_t numTrees() const override {
        return junction_trees.size();
    }

    junction_tree_map& getJunctionTrees() override {
        return junction_trees;
    }

    JunctionTreeNode& createOrGetTree(Kmer const& kmer) override;

    template<typename Archive>
    void serialize(Archive& ar) {
        ar(junction_trees, color);
    }

private:
    /// A map which stores the junction trees associated with each k-mer
    junction_tree_map junction_trees;
};

void define_MemLinkDB(py::module& m);

}

#endif //PYFROST_MEMLINKDB_H
