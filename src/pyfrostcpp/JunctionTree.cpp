#include "JunctionTree.h"

#include <utility>

namespace pyfrost {

JunctionTreeNode::JunctionTreeNode() : parent_edge(0), parent(nullptr), count(1)
{
}

JunctionTreeNode::JunctionTreeNode(char _parent_edge, JunctionTreeNode::parent_ptr_t _parent)
    : parent_edge(_parent_edge), parent(_parent), count(0)
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

    children[edge]->increment();

    return *children[edge];
}

void JunctionTreeNode::increment()
{
    ++count;
}

size_t JunctionTreeNode::getCount() const
{
    return count;
}

bool JunctionTreeNode::isLeaf() const
{
    if(children.empty()) {
        return true;
    }

    size_t child_sum = 0;
    for(auto const& it : children) {
        child_sum += it.second->getCount();
    }

    return child_sum < count;
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

    return string(output.begin(), output.end());
}

vector<size_t> JunctionTreeNode::getCoverages() {
    deque<size_t> output;
    JunctionTreeNode* curr = this;

    while(curr->parent_edge) {
        output.push_front(curr->count);
        if(curr->parent != nullptr) {
            curr = curr->parent;
        } else {
            break;
        }
    }

    return {output.begin(), output.end()};
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
            return py::make_key_iterator<py::return_value_policy::reference>(
                self.getChildren().begin(), self.getChildren().end());
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
            string choices = self.getJunctionChoices();
            return hash<string>()(choices);
        })

        .def("__eq__", [] (JunctionTreeNode& self, JunctionTreeNode& other) {
            return self.getJunctionChoices() == other.getJunctionChoices();
        })

        .def("junction_choices", &JunctionTreeNode::getJunctionChoices)
        .def("coverages", [] (JunctionTreeNode& self) {
            return as_pyarray<vector<size_t>>(move(self.getCoverages()));
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
}

}
