#ifndef PYFROST_JUNCTIONTREE_H
#define PYFROST_JUNCTIONTREE_H

#include "pyfrost.h"
#include <robin_hood.h>

#include "Kmer.h"
#include "Serialize.h"

namespace pyfrost {

class JunctionTreeNode {
public:
    using junction_node_ptr_t = std::unique_ptr<JunctionTreeNode>;
    using parent_ptr_t = JunctionTreeNode*;
    using children_t = std::map<char, junction_node_ptr_t>;

    JunctionTreeNode();
    JunctionTreeNode(char _parent_edge, parent_ptr_t parent);
    JunctionTreeNode(JunctionTreeNode&& o) = default;

    parent_ptr_t getParent();
    char getParentEdge() const;
    children_t& getChildren();
    JunctionTreeNode& addEdge(char edge);

    bool isLeaf() const;

    void prune(size_t threshold);
    size_t getCount() const;

    string getJunctionChoices();
    vector<size_t> getCoverages();

    template<typename Archive>
    void serialize(Archive& ar) {
        ar(parent_edge, count, children);
    }

    // To be called after loading from a file.
    void fixParents() {
        for(auto& child : children) {
            child.second->parent = this;
            child.second->fixParents();
        }
    }

private:
    void increment();

    char parent_edge;
    parent_ptr_t parent;
    children_t children;
    size_t count;
};

void define_JunctionTreeNode(py::module& m);

}

#endif //PYFROST_JUNCTIONTREE_H
