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

    class keys_iterator : public children_t::iterator {
    public:
        using value_type = char;

        keys_iterator() : children_t::iterator() { };
        keys_iterator(children_t::iterator o) : children_t::iterator(o) { };

        char operator*() {
            return children_t::iterator::operator*().first;
        }

        char const* operator->() {
            return &(children_t::iterator::operator->()->first);
        }
    };

    class values_iterator : public children_t::iterator {
    public:
        using value_type = JunctionTreeNode*;

        values_iterator() : children_t::iterator() { };
        values_iterator(children_t::iterator o) : children_t::iterator(o) { };

        JunctionTreeNode* operator*() {
            return children_t::iterator::operator*().second.get();
        }

        JunctionTreeNode* operator->() {
            return children_t::iterator::operator->()->second.get();
        }
    };

    JunctionTreeNode();
    JunctionTreeNode(char _parent_edge, parent_ptr_t parent);
    JunctionTreeNode(JunctionTreeNode&& o) = default;

    parent_ptr_t getParent();
    char getParentEdge() const;
    children_t& getChildren();
    JunctionTreeNode& addEdge(char edge);

    keys_iterator keys_begin() {
        return {children.begin()};
    }

    keys_iterator keys_end() {
        return {children.end()};
    }

    values_iterator values_begin() {
        return {children.begin()};
    }

    values_iterator values_end() {
        return {children.end()};
    }

    bool isLeaf() const;

    void prune(size_t threshold);
    uint16_t getCount() const;

    string getJunctionChoices();
    vector<uint16_t> getCoverages();

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

    parent_ptr_t parent;
    children_t children;
    uint16_t count;
    char parent_edge;
};

void define_JunctionTreeNode(py::module& m);

}

#endif //PYFROST_JUNCTIONTREE_H
