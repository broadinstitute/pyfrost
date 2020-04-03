#include <tuple>

#include <pybind11/pybind11.h>
#include <CompactedDBG.hpp>

#include "Pyfrost.h"

#ifndef PYFROST_NODEVIEW_H
#define PYFROST_NODEVIEW_H

namespace py = pybind11;

namespace pyfrost {

using std::pair;

class NodeView {
public:
    explicit NodeView(PyfrostCCDBG& dbg) : dbg(dbg) { }

    /**
     * This function searches for a unitig which starts with the given kmer.
     *
     * This function differs from `findKmer` because it doesn't consider internal k-mers.
     *
     * @param kmer
     */
    inline PyfrostColoredUMap findNode(Kmer const& kmer) {
        auto unitig = dbg.find(kmer, true);

        if(unitig.isEmpty) {
            throw std::out_of_range("Node not found.");
        }

        return unitig.mappingToFullUnitig();
    }

    inline PyfrostColoredUMap findNode(char const* kmer) {
        return findNode(Kmer(kmer));
    }

    /**
     * Return itself when calling with an existing kmer on unitig
     */
    inline PyfrostColoredUMap findNode(PyfrostColoredUMap const& unitig) {
        if(!unitig.isFullMapping()) {
            throw std::out_of_range("Node not found.");
        }

        return unitig.mappingToFullUnitig();
    }

    bool contains(PyfrostColoredUMap const& unitig) {
        return unitig.isFullMapping();
    }

    bool contains(Kmer const& kmer) {
        auto unitig = dbg.find(kmer, true);

        return !unitig.isEmpty;
    }

    bool contains(char const* kmer) {
        return contains(Kmer(kmer));
    }

    /**
     * Iterate over unitigs (nodes)
     */
    inline PyfrostCCDBG::iterator begin() const {
        return dbg.begin();
    }

    /**
     * End of the unitig (node) iterator
     */
    inline PyfrostCCDBG::iterator end() const {
        return dbg.end();
    }

    size_t numNodes() const {
        return dbg.size();
    }

private:
    PyfrostCCDBG& dbg;
};


/**
 * This class wraps a UnitigIterator, and adds a node data at dereference.
 *
 * There are multiple options available for iteration:
 *
 * - Return just the node: node
 * - Return a tuple with the full data dictionary: (node, full_data_dict)
 * - Return a tuple with a specific key from the data dict: (node, value_for_given_key)
 *
 * If using the last option, and the given key doesn't exist in the data dictionary, the default value can be set too
 * (defaults to None).
 */
class NodeWithDataIterator {
private:
    PyfrostCCDBG::iterator wrapped;
    bool data = false;
    py::object data_key;
    py::object default_value;

public:
    using iterator_category = std::input_iterator_tag;
    using value_type = py::object;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    explicit NodeWithDataIterator(PyfrostCCDBG::iterator const& to_wrap)
        : wrapped(to_wrap), default_value(py::none()) { }

    NodeWithDataIterator(PyfrostCCDBG::iterator const& to_wrap, bool data)
        : wrapped(to_wrap), data(data), data_key(py::none()), default_value(py::none()) { }

    NodeWithDataIterator(PyfrostCCDBG::iterator const& to_wrap, bool data, py::object data_key)
        : wrapped(to_wrap), data(data), data_key(std::move(data_key)) { }

    NodeWithDataIterator(PyfrostCCDBG::iterator const& to_wrap, bool data, py::object data_key,
                         py::object default_value)
        : wrapped(to_wrap), data(data), data_key(std::move(data_key)), default_value(std::move(default_value)) { }

    NodeWithDataIterator(NodeWithDataIterator const& o) = default;
    NodeWithDataIterator(NodeWithDataIterator&& o) = default;

    value_type operator*() {
        if(data) {
            // TODO: Populate with actual node data
            py::dict node_dict;
            if(!data_key.is_none()) {
                if(node_dict.contains(data_key)) {
                    return py::make_tuple(PyfrostColoredUMap(*wrapped), node_dict[data_key]);
                } else {
                    return py::make_tuple(PyfrostColoredUMap(*wrapped), default_value);
                }
            } else {
                return py::make_tuple(PyfrostColoredUMap(*wrapped), node_dict);
            }
        } else {
            return py::cast(PyfrostColoredUMap(*wrapped));
        }
    }

    bool operator==(NodeWithDataIterator const& o) {
        return wrapped == o.wrapped && data == o.data;
    }

    bool operator!=(NodeWithDataIterator const& o) {
        return !operator==(o);
    }

    NodeWithDataIterator& operator++() {
        ++wrapped;
        return *this;
    }

    NodeWithDataIterator operator++(int) {
        NodeWithDataIterator tmp(*this);
        ++wrapped;

        return tmp;
    }
};


class NodeDataView {
public:
    explicit NodeDataView(NodeView& node_view)
        : node_view(node_view), data_key(py::none()), default_value(py::none()) { }
    NodeDataView(NodeView& node_view, bool data)
        : node_view(node_view), data(data), data_key(py::none()), default_value(py::none()) { }
    NodeDataView(NodeView& node_view, bool data, py::object data_key)
        : node_view(node_view), data(data), data_key(std::move(data_key)), default_value(py::none()) { }
    NodeDataView(NodeView& node_view, bool data, py::object data_key, py::object default_value)
        : node_view(node_view), data(data), data_key(std::move(data_key)), default_value(std::move(default_value)) { }

    NodeDataView(NodeDataView const& o) = delete;
    NodeDataView(NodeDataView&& o) = delete;

    NodeDataView& operator=(NodeDataView& o) = delete;

    NodeWithDataIterator begin() const {
        return NodeWithDataIterator(node_view.begin(), data, data_key, default_value);
    }

    NodeWithDataIterator end() const {
        return NodeWithDataIterator(node_view.end(), data, data_key, default_value);
    }

    size_t numNodes() const {
        return node_view.numNodes();
    }

    bool contains(PyfrostColoredUMap const& unitig) {
        return node_view.contains(unitig);
    }

    bool contains(Kmer const& kmer) {
        return node_view.contains(kmer);
    }

    bool contains(char const* kmer) {
        return node_view.contains(kmer);
    }

private:
    NodeView& node_view;
    bool data = false;
    py::object data_key;
    py::object default_value;
};


void define_NodeView(py::module& m) {
    py::class_<NodeDataView>(m, "NodeDataView")
        .def("__contains__", py::overload_cast<char const*>(&NodeDataView::contains), py::is_operator())
        .def("__contains__", py::overload_cast<Kmer const&>(&NodeDataView::contains), py::is_operator())
        .def("__contains__", py::overload_cast<PyfrostColoredUMap const&>(&NodeDataView::contains), py::is_operator())

        .def("__len__", &NodeDataView::numNodes, py::is_operator())

        // This iterator already should have copies, so no need for return_value_policy::copy
        .def("__iter__", [] (NodeDataView const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    py::class_<NodeView>(m, "NodeView")
        // Access nodes with [] operator overloading
        .def("__getitem__", py::overload_cast<char const*>(&NodeView::findNode), py::is_operator())
        .def("__getitem__", py::overload_cast<Kmer const&>(&NodeView::findNode), py::is_operator())
        .def("__getitem__", py::overload_cast<PyfrostColoredUMap const&>(&NodeView::findNode), py::is_operator())

        .def("__contains__", py::overload_cast<char const*>(&NodeView::contains), py::is_operator())
        .def("__contains__", py::overload_cast<Kmer const&>(&NodeView::contains), py::is_operator())
        .def("__contains__", py::overload_cast<PyfrostColoredUMap const&>(&NodeView::contains), py::is_operator())

        .def("__len__", &NodeView::numNodes, py::is_operator())

        // Call without parameters just returns the node iterator
        .def("__call__", [] (NodeView const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.begin(), self.end());
        })

        // With data arguments we return a NodeDataView
        .def("__call__", [] (NodeView& self, py::object const& data) {
            try {
                bool data_bool = py::cast<bool>(data);

                return new NodeDataView(self, data_bool);
            } catch(py::cast_error const&) {
                return new NodeDataView(self, true, data);
            }
        }, py::keep_alive<0, 1>())
        .def("__call__", [] (NodeView& self, py::kwargs const& kwargs) {
            py::object default_value = py::none();
            if(kwargs.contains("default")) {
                default_value = kwargs["default"];
            }

            py::object data = py::none();
            bool data_bool = false;
            if(kwargs.contains("data")) {
                try {
                    data_bool = py::cast<bool>(kwargs["data"]);

                    return new NodeDataView(self, data_bool);
                } catch(py::cast_error const&) {
                    return new NodeDataView(self, true, kwargs["data"], default_value);
                }
            }

            return new NodeDataView(self);
        }, py::keep_alive<0, 1>())

        // Iterate over all nodes
        .def("__iter__", [](NodeView const& self) {
            return py::make_iterator<py::return_value_policy::copy>(self.begin(), self.end());
        }, py::keep_alive<0, 1>());
}


}

#endif //PYFROST_NODEVIEW_H
