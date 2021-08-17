#ifndef PYFROST_JUNCTIONTREE_H
#define PYFROST_JUNCTIONTREE_H

#include "Kmer.h"
#include <robin_hood.h>

namespace pyfrost {

class JunctionTreeNode : public std::enable_shared_from_this<JunctionTreeNode> {
public:
    using junction_node_ptr_t = std::shared_ptr<JunctionTreeNode>;
    using parent_ptr_t = std::weak_ptr<JunctionTreeNode>;
    using children_t = std::map<char, junction_node_ptr_t>;

    JunctionTreeNode();
    JunctionTreeNode(char _parent_edge, parent_ptr_t parent);
    JunctionTreeNode(JunctionTreeNode const& o) = default;
    JunctionTreeNode(JunctionTreeNode&& o) = default;

    parent_ptr_t getParent();
    char getParentEdge() const;
    children_t& getChildren();
    junction_node_ptr_t addEdge(char edge);

    bool isLeaf() const;

    void prune(size_t threshold);
    size_t getCount() const;

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
