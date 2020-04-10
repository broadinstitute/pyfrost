#include "Pyfrost.h"
#include "UnitigDataProxy.h"

#ifndef PYFROST_NODEVIEW_H
#define PYFROST_NODEVIEW_H

namespace py = pybind11;

namespace pyfrost {

class NodeDataView;

class NodeView {
public:
    friend class NodeDataView;

    explicit NodeView(PyfrostCCDBG& dbg) : dbg(dbg) { }

    /**
     * This function searches for a unitig which starts with the given kmer.
     *
     * This function differs from `findKmer` because it doesn't consider internal k-mers.
     *
     * @param kmer
     */
    inline UnitigDataProxy findNode(Kmer const& kmer) {
        auto unitig = dbg.find(kmer, true);

        if(unitig.isEmpty) {
            throw std::out_of_range("Node not found.");
        }

        unitig = unitig.mappingToFullUnitig();
        return UnitigDataProxy(unitig);
    }

    inline UnitigDataProxy findNode(char const* kmer) {
        return findNode(Kmer(kmer));
    }

    /**
     * Return itself when calling with an existing kmer on unitig
     */
    inline UnitigDataProxy findNode(PyfrostColoredUMap const& unitig) {
        if(!unitig.isFullMapping()) {
            throw std::out_of_range("Node not found.");
        }

        return UnitigDataProxy(unitig);
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
            auto data_proxy = UnitigDataProxy(*wrapped);

            if(!data_key.is_none()) {
                auto key = data_key.cast<std::string>();
                if(data_proxy.contains(key)) {
                    return py::make_tuple(PyfrostColoredUMap(*wrapped), data_proxy.getData(key));
                } else {
                    return py::make_tuple(PyfrostColoredUMap(*wrapped), default_value);
                }
            } else {
                return py::make_tuple(PyfrostColoredUMap(*wrapped), data_proxy);
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

    py::object get(PyfrostColoredUMap const& unitig) {
        if(data) {
            auto data_proxy = UnitigDataProxy(unitig);

            if(!data_key.is_none()) {
                auto key = data_key.cast<std::string>();
                if(data_proxy.contains(key)) {
                    return data_proxy.getData(key);
                } else {
                    return default_value;
                }
            } else {
                return py::cast(data_proxy);
            }
        } else {
            return py::cast(unitig);
        }
    }

    inline py::object get(Kmer const& kmer) {
        auto unitig = node_view.dbg.find(kmer);
        return get(unitig);
    }

    inline py::object get(char const* kmer) {
        return get(Kmer(kmer));
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


void define_NodeView(py::module& m);

}

#endif //PYFROST_NODEVIEW_H
