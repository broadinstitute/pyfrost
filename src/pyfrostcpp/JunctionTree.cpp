#include "JunctionTree.h"

#include <utility>

namespace pyfrost {

JunctionTreeNode::JunctionTreeNode() : parent(nullptr), parent_edge(0)
{
}

JunctionTreeNode::JunctionTreeNode(char _parent_edge, JunctionTreeNode::parent_ptr_t _parent)
    : parent(_parent), parent_edge(_parent_edge)
{
}

JunctionTreeNode::parent_ptr_t JunctionTreeNode::getParent()
{
    return parent;
}

char JunctionTreeNode::getParentEdge() const {
    return parent_edge;
}

JunctionTreeNode::children_t& JunctionTreeNode::getChildren()
{
    return children;
}


JunctionTreeNode& JunctionTreeNode::addEdge(char edge)
{
    if(!(edge == 'A' || edge == 'C' || edge =='G' || edge == 'T')) {
        throw std::runtime_error("Invalid edge, must be one of A, C, T or G");
    }

    auto it = children.find(edge);
    if(it == children.end()) {
        children.emplace(edge, std::make_unique<JunctionTreeNode>(edge, this));
    }

    return *children[edge];
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

    while(curr->parent_edge) {
        output.push_front(curr->parent_edge);
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

    while(curr->parent_edge) {
        output.push_front(curr->getCount());
        if(curr->parent != nullptr) {
            curr = curr->parent;
        } else {
            break;
        }
    }

    return {output.begin(), output.end()};
}

// JunctionTreeNodeWithCov
uint16_t JunctionTreeNodeWithCov::getCount() const
{
    return count;
}

void JunctionTreeNodeWithCov::increment()
{
    ++count;
}

JunctionTreeNode& JunctionTreeNodeWithCov::addEdge(char edge)
{
    if(!(edge == 'A' || edge == 'C' || edge =='G' || edge == 'T')) {
        throw std::runtime_error("Invalid edge, must be one of A, C, T or G");
    }

    auto it = children.find(edge);
    if(it == getChildren().end()) {
        children.emplace(edge, std::make_unique<JunctionTreeNodeWithCov>(edge, this));
    }

    dynamic_cast<JunctionTreeNodeWithCov*>(children[edge].get())->increment();

    return *children[edge];
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

            char parent = self.getParentEdge();
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

        .def("prune", &JunctionTreeNode::prune)
        .def("is_leaf", &JunctionTreeNode::isLeaf)

        .def_property_readonly("count", &JunctionTreeNode::getCount)
        .def_property_readonly("parent", [] (JunctionTreeNode& self) -> py::object {
            if(self.getParent() != nullptr) {
                return py::cast(self.getParent());
            }

            return py::none();
        }, py::return_value_policy::reference)
        .def_property_readonly("parent_edge", [] (JunctionTreeNode& self) -> py::object {
            if(self.getParentEdge() != 0) {
                return py::cast(self.getParentEdge());
            }

            return py::none();
        }, py::return_value_policy::reference);

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    py_JunctionTreeNode.attr("__bases__") = py::make_tuple(Mapping).attr("__add__")(py_JunctionTreeNode.attr("__bases__"));

    auto py_JunctionTreeNodeWithCov = py::class_<JunctionTreeNodeWithCov, JunctionTreeNode>(
        m, "JunctionTreeNodeWithCov")
        .def("is_leaf", &JunctionTreeNodeWithCov::isLeaf)
        .def_property_readonly("count", &JunctionTreeNodeWithCov::getCount);

}

}
