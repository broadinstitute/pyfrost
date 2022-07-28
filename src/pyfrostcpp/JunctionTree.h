#ifndef PYFROST_JUNCTIONTREE_H
#define PYFROST_JUNCTIONTREE_H

#include "pyfrost.h"
#include <robin_hood.h>

#include "Kmer.h"
#include "Serialize.h"

#include <cereal/types/base_class.hpp>
#include <cereal/archives/binary.hpp>

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
    virtual JunctionTreeNode& addEdge(char edge);

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

    virtual bool isLeaf() const;

    void prune(size_t threshold);
    virtual uint16_t getCount() const;

    string getJunctionChoices();
    vector<uint16_t> getCoverages();

    template<typename Archive>
    void serialize(Archive& ar) {
        ar(parent_edge, children);
    }

    // To be called after loading from a file.
    void fixParents() {
        for(auto& child : children) {
            child.second->parent = this;
            child.second->fixParents();
        }
    }

protected:
    parent_ptr_t parent;
    children_t children;
    char parent_edge;
};


class JunctionTreeNodeWithCov : public JunctionTreeNode {
public:
    JunctionTreeNodeWithCov() : count(1), JunctionTreeNode() { }
    JunctionTreeNodeWithCov(char _parent_edge, parent_ptr_t parent) :
        count(0), JunctionTreeNode(_parent_edge, parent) { }
    JunctionTreeNodeWithCov(JunctionTreeNodeWithCov&& o) = default;

    JunctionTreeNode& addEdge(char edge) override;
    uint16_t getCount() const override;

    template<typename Archive>
    void serialize(Archive& ar) {
        ar(cereal::base_class<JunctionTreeNode>(this), count);
    }

private:
    void increment();

    uint16_t count;
};

void define_JunctionTreeNode(py::module& m);

}

CEREAL_REGISTER_TYPE(pyfrost::JunctionTreeNodeWithCov)
CEREAL_REGISTER_POLYMORPHIC_RELATION(pyfrost::JunctionTreeNode, pyfrost::JunctionTreeNodeWithCov)

#endif //PYFROST_JUNCTIONTREE_H
