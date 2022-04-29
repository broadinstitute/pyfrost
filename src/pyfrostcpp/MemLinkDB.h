#ifndef PYFROST_MEMLINKDB_H
#define PYFROST_MEMLINKDB_H

#include "pyfrost.h"
#include "JunctionTree.h"
#include "LinkDB.h"
#include "Serialize.h"

#include <cereal/archives/binary.hpp>

namespace pyfrost {

/**
 * Store all links in memory.
 *
 * @tparam T The type of JunctionTreeNode, can also be JunctionTreeNodeWithCov
 */
template<typename T>
class MemLinkDB : public LinkDB {
public:
    MemLinkDB() : LinkDB() { }
    explicit MemLinkDB(size_t _color) : LinkDB(_color) { }
    MemLinkDB(MemLinkDB&& o) noexcept = default;

    MemLinkDB& operator=(MemLinkDB const& o) = default;

    bool hasLinks(Kmer const& kmer) const override {
        return junction_trees.contains(kmer);
    }

    T& getLinks(Kmer const& kmer) override;

    size_t numTrees() const override {
        return junction_trees.size();
    }

    junction_tree_map& getJunctionTrees() override {
        return junction_trees;
    }

    T& createOrGetTree(Kmer const& kmer) override;

    template<typename Archive>
    void serialize(Archive& ar) {
        ar(junction_trees, color);
    }

private:
    /// A map which stores the junction trees associated with each k-mer
    junction_tree_map junction_trees;
};


template<typename T>
T& MemLinkDB<T>::createOrGetTree(Kmer const &kmer) {
    auto it = junction_trees.find(kmer);
    if (it == junction_trees.end()) {
        junction_trees.emplace(kmer, std::make_unique<T>());
        return *dynamic_cast<T*>(junction_trees[kmer].get());
    } else {
        return *dynamic_cast<T*>(it->second.get());
    }
}

template<typename T>
T& MemLinkDB<T>::getLinks(const Kmer &kmer) {
    return *dynamic_cast<T*>(junction_trees.at(kmer).get());
}

template<typename T>
void define_MemLinkDB(py::module& m, const char* name) {
    auto py_MemLinkDB = py::class_<MemLinkDB<T>, LinkDB>(m, name)
        .def(py::init())
        .def(py::init<size_t>())
        .def("get_links", [] (MemLinkDB<T>& self, Kmer const& kmer) -> JunctionTreeNode& {
            return self.getLinks(kmer);
        }, py::return_value_policy::reference)

        .def("__contains__", &MemLinkDB<T>::hasLinks)

        .def("__getitem__", [] (MemLinkDB<T>& self, const Kmer& kmer) -> JunctionTreeNode& {
            return self.getLinks(kmer);
        }, py::return_value_policy::reference)

        .def("__len__", [] (MemLinkDB<T>& self) {
            return self.numTrees();
        })
        .def("__iter__", [] (MemLinkDB<T>& self) {
            return py::make_key_iterator(self.getJunctionTrees().begin(), self.getJunctionTrees().end());
        }, py::keep_alive<0, 1>())

        .def("save", [] (MemLinkDB<T>& self, string const& filepath) {
            std::ofstream ofile(filepath);
            cereal::BinaryOutputArchive archive(ofile);

            archive(cereal::make_nvp("memlinkdb", self));
        })
        .def_static("from_file", [] (string const& filepath) {
            MemLinkDB<T> linkdb;

            std::ifstream ifile(filepath);
            cereal::BinaryInputArchive archive(ifile);

            archive(linkdb);
            linkdb.fixTreeParents();

            return linkdb;
        });

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    auto Set = py::module::import("collections.abc").attr("Set");
    py_MemLinkDB.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(py_MemLinkDB.attr("__bases__"));
}

}

#endif //PYFROST_MEMLINKDB_H
