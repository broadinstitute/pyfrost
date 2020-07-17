
#ifndef PYFROST_NODEVIEW_H
#define PYFROST_NODEVIEW_H

#include "Pyfrost.h"
#include "UnitigDataProxy.h"
#include "NodeIterator.h"

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
     * @param kmer
     */
    inline UnitigDataProxy findNode(Kmer const& kmer) {
        auto unitig = dbg.find(kmer, true).mappingToFullUnitig();

        if(unitig.isEmpty) {
            throw std::out_of_range("Node not found.");
        }

        return UnitigDataProxy(unitig);
    }

    inline UnitigDataProxy findNode(char const* kmer) {
        return findNode(Kmer(kmer));
    }

    /**
     * Return itself when calling with an existing kmer on unitig
     */
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
    inline auto begin() const {
        return NodeIterator<PyfrostCCDBG::iterator>(&dbg, dbg.begin(), false);
    }

    /**
     * End of the unitig (node) iterator
     */
    inline auto end() const {
        return NodeIterator<PyfrostCCDBG::iterator>(&dbg, dbg.end(), false);
    }

    size_t numNodes() const {
        return dbg.size();
    }

private:
    PyfrostCCDBG& dbg;
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

    auto begin() const {
        return NodeWithDataIterator<PyfrostCCDBG::iterator>(node_view.begin(), data, data_key, default_value);
    }

    auto end() const {
        return NodeWithDataIterator<PyfrostCCDBG::iterator>(node_view.end(), data, data_key, default_value);
    }

    size_t numNodes() const {
        return node_view.numNodes();
    }

    inline py::object get(Kmer const& kmer) {
        auto unitig = node_view.dbg.find(kmer, true).mappingToFullUnitig();

        if(unitig.isEmpty) {
            throw py::key_error("Unitig does not exists");
        }

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

    inline py::object get(char const* kmer) {
        return get(Kmer(kmer));
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
