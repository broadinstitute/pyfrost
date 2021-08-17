#include "MemLinkDB.h"

namespace pyfrost {

std::shared_ptr<JunctionTreeNode> MemLinkDB::createOrGetTree(Kmer const &kmer) {
    auto it = junction_trees.find(kmer);
    if (it == junction_trees.end()) {
        junction_trees.emplace(kmer, std::make_shared<JunctionTreeNode>());
        return junction_trees[kmer];
    } else {
        return it->second;
    }
}

std::shared_ptr<JunctionTreeNode> MemLinkDB::getLinks(const Kmer &kmer) {
    return junction_trees.at(kmer);
}

void define_MemLinkDB(py::module& m) {
    auto py_MemLinkDB = py::class_<MemLinkDB, LinkDB>(m, "MemLinkDB")
        .def(py::init())
        .def("get_links", [] (MemLinkDB& self, Kmer const& kmer) {
            return self.getLinks(kmer);
        }, py::keep_alive<0, 1>())

        .def("__contains__", &MemLinkDB::hasLinks)

        .def("__getitem__", [] (MemLinkDB& self, const Kmer& kmer) {
            return self.getLinks(kmer);
        }, py::keep_alive<0, 1>())

        .def("__len__", [] (MemLinkDB& self) {
            return self.numTrees();
        })
        .def("__iter__", [] (MemLinkDB& self) {
            return py::make_key_iterator(self.getJunctionTrees().begin(), self.getJunctionTrees().end());
        }, py::keep_alive<0, 1>());

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    auto Set = py::module::import("collections.abc").attr("Set");
    py_MemLinkDB.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(py_MemLinkDB.attr("__bases__"));
}

}
