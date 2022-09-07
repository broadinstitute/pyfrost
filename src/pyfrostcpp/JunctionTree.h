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
    JunctionTreeNode(char _label, parent_ptr_t parent);
    JunctionTreeNode(JunctionTreeNode&& o) = default;

    parent_ptr_t getParent();
    char getLabel() const;
    children_t& getChildren();
    virtual JunctionTreeNode& addOrGetChild(char child);
    virtual JunctionTreeNode& addOrIncrementChild(char child);

    // No coverage stored in this class, but here to keep API consistent
    virtual void increment() { }
    virtual void increment(uint32_t num) { }

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

    virtual void merge(JunctionTreeNode& other);

    template<typename Archive>
    void serialize(Archive& ar) {
        ar(label, children);
    }

    // To be called after loading from a file.
    void fixParents() {
        for(auto& child : children) {
            child.second->parent = this;
            child.second->fixParents();
        }
    }

public:

    parent_ptr_t parent;
    children_t children;
    char label;
};


class JunctionTreeNodeWithCov : public JunctionTreeNode {
public:
    JunctionTreeNodeWithCov() : count(1), JunctionTreeNode() { }
    JunctionTreeNodeWithCov(char label, parent_ptr_t parent) :
        count(0), JunctionTreeNode(label, parent) { }
    JunctionTreeNodeWithCov(JunctionTreeNodeWithCov&& o) = default;

    JunctionTreeNode& addOrGetChild(char label) override;
    JunctionTreeNode& addOrIncrementChild(char label) override;
    uint16_t getCount() const override;

    void merge(JunctionTreeNode& other) override;

    void increment() override;
    void increment(uint32_t num) override;

    template<typename Archive>
    void serialize(Archive& ar) {
        ar(cereal::base_class<JunctionTreeNode>(this), count);
    }

private:
    uint16_t count;
};

void define_JunctionTreeNode(py::module& m);

}

CEREAL_REGISTER_TYPE(pyfrost::JunctionTreeNodeWithCov)
CEREAL_REGISTER_POLYMORPHIC_RELATION(pyfrost::JunctionTreeNode, pyfrost::JunctionTreeNodeWithCov)

#endif //PYFROST_JUNCTIONTREE_H
