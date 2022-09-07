#include "JunctionTree.h"

#include <utility>
#include <queue>

namespace pyfrost {

JunctionTreeNode::JunctionTreeNode() : parent(nullptr), label(0)
{
}

JunctionTreeNode::JunctionTreeNode(char _label, JunctionTreeNode::parent_ptr_t _parent)
    : parent(_parent), label(_label)
{
}

JunctionTreeNode::parent_ptr_t JunctionTreeNode::getParent()
{
    return parent;
}

char JunctionTreeNode::getLabel() const {
    return label;
}

JunctionTreeNode::children_t& JunctionTreeNode::getChildren()
{
    return children;
}

JunctionTreeNode& JunctionTreeNode::addOrGetChild(char child)
{
    if(!(child == 'A' || child == 'C' || child == 'G' || child == 'T')) {
        throw std::runtime_error("Invalid child, must be one of A, C, T or G");
    }

    auto it = children.find(child);
    if(it == getChildren().end()) {
        children.emplace(child, std::make_unique<JunctionTreeNode>(child, this));
    }

    return *children[child];
}

/**
 * Because this JunctionTreeNode doesn't store coverage, this is just a wrapper around addOrGetChild
 */
JunctionTreeNode& JunctionTreeNode::addOrIncrementChild(char child)
{
    return addOrGetChild(child);
}


uint16_t JunctionTreeNode::getCount() const
{
    return 1;
}

bool JunctionTreeNode::isLeaf() const
{
    return children.empty();
}

void JunctionTreeNode::prune(size_t threshold) {
    std::vector<char> to_delete;
    for(auto& pair : children) {
        if(pair.second->getCount() < threshold) {
            to_delete.push_back(pair.first);
        }
    }

    for(char edge : to_delete) {
        children.erase(edge);
    }

    for(auto& pair : children) {
        pair.second->prune(threshold);
    }
}

string JunctionTreeNode::getJunctionChoices() {
    deque<char> output;
    JunctionTreeNode* curr = this;

    while(curr->label) {
        output.push_front(curr->label);
        if(curr->parent != nullptr) {
            curr = curr->parent;
        } else {
            break;
        }
    }

    return {output.begin(), output.end()};
}

vector<uint16_t> JunctionTreeNode::getCoverages() {
    deque<size_t> output;
    JunctionTreeNode* curr = this;

    while(curr->label) {
        output.push_front(curr->getCount());
        if(curr->parent != nullptr) {
            curr = curr->parent;
        } else {
            break;
        }
    }

    return {output.begin(), output.end()};
}


void JunctionTreeNode::merge(JunctionTreeNode& other)
{
    std::queue<pair<JunctionTreeNode*, JunctionTreeNode*>> q;
    q.push(make_pair(this, &other));

    while(!q.empty()) {
        JunctionTreeNode* target = q.front().first;
        JunctionTreeNode* source = q.front().second;
        q.pop();

        for(auto& other_child : source->getChildren()) {
            auto& target_child = target->addOrGetChild(other_child.first);
            q.push(make_pair(&target_child, other_child.second.get()));
        }
    }
}

// JunctionTreeNodeWithCov
uint16_t JunctionTreeNodeWithCov::getCount() const
{
    return count;
}

void JunctionTreeNodeWithCov::increment()
{
    if(count < numeric_limits<uint16_t>::max()) {  // Don't overflow
        ++count;
    }
}

void JunctionTreeNodeWithCov::increment(uint32_t num)
{
    count = min(static_cast<uint32_t>(numeric_limits<uint16_t>::max()), num + count);
}

JunctionTreeNode& JunctionTreeNodeWithCov::addOrGetChild(char label)
{
    if(!(label == 'A' || label == 'C' || label =='G' || label == 'T')) {
        throw std::runtime_error("Invalid label, must be one of A, C, T or G");
    }

    auto it = children.find(label);
    if(it == getChildren().end()) {
        children.emplace(label, std::make_unique<JunctionTreeNodeWithCov>(label, this));
    }

    return *children[label];
}

JunctionTreeNode& JunctionTreeNodeWithCov::addOrIncrementChild(char label)
{
    auto& child = addOrGetChild(label);
    dynamic_cast<JunctionTreeNodeWithCov*>(&child)->increment();

    return child;
}

void JunctionTreeNodeWithCov::merge(JunctionTreeNode& other)
{
    std::queue<pair<JunctionTreeNode*, JunctionTreeNode*>> q;
    q.push(make_pair(this, &other));

    while(!q.empty()) {
        JunctionTreeNode* target = q.front().first;
        JunctionTreeNode* source = q.front().second;
        q.pop();

        target->increment(source->getCount());

        for(auto& other_child : source->getChildren()) {
            auto& target_child = target->addOrGetChild(other_child.first);
            q.push(make_pair(&target_child, other_child.second.get()));
        }
    }
}


void define_JunctionTreeNode(py::module& m) {
    auto py_JunctionTreeNode = py::class_<JunctionTreeNode>(m, "JunctionTreeNode")
        .def("__getitem__", [] (JunctionTreeNode& self, char edge) {
            return self.getChildren().at(edge).get();
        }, py::return_value_policy::reference)
        .def("__len__", [] (JunctionTreeNode& self) {
            return self.getChildren().size();
        })
        .def("__iter__", [] (JunctionTreeNode& self) {
            return py::make_iterator(self.keys_begin(), self.keys_end());
        }, py::keep_alive<0, 1>())

        .def("keys", [] (JunctionTreeNode& self) {
            return py::make_iterator(self.keys_begin(), self.keys_end());
        }, py::keep_alive<0, 1>())

        .def("values", [] (JunctionTreeNode& self) {
            return py::make_iterator<py::return_value_policy::reference>(self.values_begin(), self.values_end());
        }, py::keep_alive<0, 1>())

        .def("__repr__", [] (JunctionTreeNode& self) {
            std::stringstream sstream;

            char parent = self.getLabel();
            parent = parent == 0 ? '-' : parent;
            sstream << "<JunctionTreeNode parent=" << parent << " count=" << self.getCount();
            sstream << " children=[";
            bool first = true;
            for(auto const& it : self.getChildren()) {
                if(!first) {
                    sstream << ", ";
                }
                sstream << it.first;
                first = false;
            }

            sstream << "]>";

            return sstream.str();
        })

        .def("__hash__", [] (JunctionTreeNode& self) {
            // Ownership of junction tree nodes is well defined in a tree, so we just check if the addresses match
            // I.e., we expect that link databases are loaded only once, and each tree node then has a fixed memory
            // address.
            return hash<JunctionTreeNode*>()(&self);
        })

        .def("__eq__", [] (JunctionTreeNode& self, JunctionTreeNode& other) {
            // Ownership of junction tree nodes is well defined in a tree, so we just check if the addresses match
            // I.e., we expect that link databases are loaded only once, and each tree node then has a fixed memory
            // address.
            return &self == &other;
        })

        .def("junction_choices", [] (JunctionTreeNode& self) {
            return py::bytes(self.getJunctionChoices());
        })
        .def("coverages", [] (JunctionTreeNode& self) {
            return as_pyarray<vector<uint16_t>>(move(self.getCoverages()));
        })

        .def("merge", &JunctionTreeNode::merge)
        .def("prune", &JunctionTreeNode::prune)
        .def("is_leaf", &JunctionTreeNode::isLeaf)

        .def_property_readonly("count", &JunctionTreeNode::getCount)
        .def_property_readonly("parent", [] (JunctionTreeNode& self) -> py::object {
            if(self.getParent() != nullptr) {
                return py::cast(self.getParent());
            }

            return py::none();
        }, py::return_value_policy::reference)
        .def_property_readonly("label", [] (JunctionTreeNode& self) -> py::object {
            if(self.getLabel() != 0) {
                return py::cast(self.getLabel());
            }

            return py::none();
        }, py::return_value_policy::reference);

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    py_JunctionTreeNode.attr("__bases__") = py::make_tuple(Mapping).attr("__add__")(py_JunctionTreeNode.attr("__bases__"));

    auto py_JunctionTreeNodeWithCov = py::class_<JunctionTreeNodeWithCov, JunctionTreeNode>(
        m, "JunctionTreeNodeWithCov")

        .def("merge", &JunctionTreeNodeWithCov::merge)

        .def("is_leaf", &JunctionTreeNodeWithCov::isLeaf)
        .def_property_readonly("count", &JunctionTreeNodeWithCov::getCount);

}

}
