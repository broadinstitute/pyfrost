#include "MemLinkDB.h"

#include <cereal/archives/binary.hpp>

namespace pyfrost {

JunctionTreeNode& MemLinkDB::createOrGetTree(Kmer const &kmer) {
    auto it = junction_trees.find(kmer);
    if (it == junction_trees.end()) {
        junction_trees.emplace(kmer, std::make_unique<JunctionTreeNode>());
        return *junction_trees[kmer];
    } else {
        return *it->second;
    }
}

JunctionTreeNode& MemLinkDB::getLinks(const Kmer &kmer) {
    return *junction_trees.at(kmer);
}

void define_MemLinkDB(py::module& m) {
    auto py_MemLinkDB = py::class_<MemLinkDB, LinkDB>(m, "MemLinkDB")
        .def(py::init())
        .def(py::init<size_t>())
        .def("get_links", [] (MemLinkDB& self, Kmer const& kmer) -> JunctionTreeNode& {
            return self.getLinks(kmer);
        }, py::return_value_policy::reference)

        .def("__contains__", &MemLinkDB::hasLinks)

        .def("__getitem__", [] (MemLinkDB& self, const Kmer& kmer) -> JunctionTreeNode& {
            return self.getLinks(kmer);
        }, py::return_value_policy::reference)

        .def("__len__", [] (MemLinkDB& self) {
            return self.numTrees();
        })
        .def("__iter__", [] (MemLinkDB& self) {
            return py::make_key_iterator(self.getJunctionTrees().begin(), self.getJunctionTrees().end());
        }, py::keep_alive<0, 1>())

        .def("save", [] (MemLinkDB& self, string const& filepath) {
            std::ofstream ofile(filepath);
            cereal::BinaryOutputArchive archive(ofile);

            archive(cereal::make_nvp("memlinkdb", self));
        })
        .def_static("from_file", [] (string const& filepath) {
            MemLinkDB linkdb;

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
