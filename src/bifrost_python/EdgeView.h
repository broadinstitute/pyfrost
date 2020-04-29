#ifndef PYFROST_EDGEVIEW_H
#define PYFROST_EDGEVIEW_H

#include <numeric>
#include <utility>

#include "Pyfrost.h"
#include "NodeBunchIter.h"

namespace pyfrost {

py::dict makeEdgeDataDict(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2);

/**
 * Iterate over all edges for a given set of nodes (as specified by the node iterator).
 *
 * By using a templated node iterator type, you can for example both use a std::vector::iterator to specify the
 * nodes, or use the PyfrostCCDBG::iterator to iterate over all nodes in the graph.
 *
 * \tparam T Node iterator type, can be used to iterate over different node containers. Dereference of this iterator
 *           type should return a UnitigMap
 */
template<typename T, bool outgoing=true>
class EdgeIterator {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = py::object;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    EdgeIterator() = default;

    EdgeIterator(T const& begin, T const& end) : node_iter(begin), node_end(end) {
        // find a node with outgoing or incoming edges
        while(node_iter != node_end) {
            size_t num = (outgoing ? node_iter->getSuccessors().cardinality() :
                          node_iter->getPredecessors().cardinality());

            if(num > 0) {
                edge_iter = outgoing ? node_iter->getSuccessors().begin() : node_iter->getPredecessors().begin();
                edge_end = outgoing ? node_iter->getSuccessors().end() : node_iter->getPredecessors().end();

                break;
            } else {
                ++node_iter;
            }
        }
    }

    EdgeIterator(EdgeIterator const& o) = default;
    EdgeIterator(EdgeIterator&& o) = default;

    EdgeIterator& operator++();
    EdgeIterator operator++(int);

    value_type operator*();

    bool operator==(EdgeIterator const& o) const;
    bool operator!=(EdgeIterator const& o) const;

private:
    T node_iter;
    T node_end;
    PyfrostColoredUMap::neighbor_iterator edge_iter;
    PyfrostColoredUMap::neighbor_iterator edge_end;
};

template<typename T, bool outgoing>
EdgeIterator<T, outgoing>& EdgeIterator<T, outgoing>::operator++() {
    if(node_iter->isEmpty || node_iter == node_end) {
        return *this;
    }

    do {
        // Move to next edge if possible
        ++edge_iter;

        if(edge_iter == edge_end) {
            // All edges of current node visited, move to next node
            ++node_iter;

            if(node_iter == node_end) {
                break;
            }

            if(outgoing) {
                edge_iter = node_iter->getSuccessors().begin();
                edge_end = node_iter->getSuccessors().end();
            } else {
                edge_iter = node_iter->getPredecessors().begin();
                edge_end = node_iter->getPredecessors().end();
            }
        }
    } while (edge_iter == edge_end);

    return *this;
}

template<typename T, bool outgoing>
EdgeIterator<T, outgoing> EdgeIterator<T, outgoing>::operator++(int) {
    EdgeIterator<T, outgoing> tmp(*this);

    operator++();
    return tmp;
}

template<typename T, bool outgoing>
typename EdgeIterator<T, outgoing>::value_type EdgeIterator<T, outgoing>::operator*() {
    if(outgoing) {
        return py::make_tuple(*node_iter, *edge_iter);
    } else {
        return py::make_tuple(*edge_iter, *node_iter);
    }
}

template<typename T, bool outgoing>
bool EdgeIterator<T, outgoing>::operator==(EdgeIterator<T, outgoing> const& o) const {
    if(node_iter == node_end && o.node_iter == o.node_end) {
        return true;
    } else if(node_iter->isEmpty && o.node_iter->isEmpty) {
        return true;
    } else {
        return node_iter == o.node_iter && edge_iter == o.edge_iter;
    }
}

template<typename T, bool outgoing>
bool EdgeIterator<T, outgoing>::operator!=(EdgeIterator<T, outgoing> const& o) const {
    return !operator==(o);
}

//
// EdgeWithDataIterator
// --------------------
//

/**
 * Wraps the EdgeIterator and adds the data dict on dereference (or a specified key).
 *
 * @tparam T
 * @tparam outgoing
 */
template<typename T, bool outgoing=true>
class EdgeWithDataIterator {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = py::object;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    EdgeWithDataIterator(EdgeIterator<T, outgoing> const& to_wrap, py::object data, py::object default_value) :
        wrapped(to_wrap), data(std::move(data)), default_value(std::move(default_value))
        { }

    EdgeWithDataIterator(EdgeWithDataIterator const& o) = default;
    EdgeWithDataIterator(EdgeWithDataIterator&& o) = default;

    value_type operator*();

    bool operator==(EdgeWithDataIterator const& o);
    bool operator!=(EdgeWithDataIterator const& o);

    EdgeWithDataIterator& operator++();
    EdgeWithDataIterator operator++(int);

private:
    EdgeIterator<T, outgoing> wrapped;
    py::object data;
    py::object default_value;

};

template<typename T, bool outgoing>
typename EdgeWithDataIterator<T, outgoing>::value_type EdgeWithDataIterator<T, outgoing>::operator*() {
    if(data.is_none() || (py::isinstance<py::bool_>(data) && !data.cast<bool>())) {
        return *wrapped;
    } else {
        auto edge = py::cast<py::tuple>(*wrapped);
        auto data_dict = makeEdgeDataDict(
            py::cast<PyfrostColoredUMap>(edge[0]),
            py::cast<PyfrostColoredUMap>(edge[1])
        );

        if(py::isinstance<py::bool_>(data)) {
            return py::make_tuple(edge[0], edge[1], data_dict);
        } else {
            auto key = data.cast<std::string>();
            if(data_dict.contains(key)) {
                return py::make_tuple(edge[0], edge[1], data_dict[key.c_str()]);
            } else {
                return py::make_tuple(edge[0], edge[1], default_value);
            }
        }
    }
}

template<typename T, bool outgoing>
bool EdgeWithDataIterator<T, outgoing>::operator==(const EdgeWithDataIterator<T, outgoing> &o) {
    return wrapped == o.wrapped && data.equal(o.data);
}

template<typename T, bool outgoing>
bool EdgeWithDataIterator<T, outgoing>::operator!=(const EdgeWithDataIterator<T, outgoing> &o) {
    return !operator==(o);
}

template<typename T, bool outgoing>
EdgeWithDataIterator<T, outgoing>& EdgeWithDataIterator<T, outgoing>::operator++() {
    ++wrapped;

    return *this;
}

template<typename T, bool outgoing>
EdgeWithDataIterator<T, outgoing> EdgeWithDataIterator<T, outgoing>::operator++(int) {
    EdgeWithDataIterator<T, outgoing> tmp(*this);

    operator++();
    return tmp;
}

template<bool outgoing> class EdgeDataView;

/**
 * Class to obtain metadata of edges or to iterate over all of them
 * @tparam outgoing Whether to access outgoing edges or the incoming edges
 */
template<bool outgoing>
class EdgeView {
public:
    friend class EdgeDataView<outgoing>;

    explicit EdgeView(PyfrostCCDBG& dbg);
    EdgeView(EdgeView const& o) = default;
    EdgeView(EdgeView&& o) = default;

    py::dict get(py::tuple const& edge);
    py::dict get(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2);

    EdgeIterator<PyfrostCCDBG::iterator, outgoing> begin() const;
    EdgeIterator<PyfrostCCDBG::iterator, outgoing> end() const;

    bool contains(py::tuple const& edge);

    virtual size_t size() const;

private:
    PyfrostColoredUMap findNode(py::object const& n);
    bool isSuccessor(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2);

    PyfrostCCDBG& dbg;
};


template<bool outgoing>
EdgeView<outgoing>::EdgeView(PyfrostCCDBG& dbg) : dbg(dbg) { }

template<bool outgoing>
PyfrostColoredUMap EdgeView<outgoing>::findNode(py::object const& n) {
    PyfrostColoredUMap node;

    if (py::isinstance<py::str>(n)) {
        auto str = n.cast<std::string>();
        node = dbg.find(Kmer(str.c_str()), true).mappingToFullUnitig();
    } else if (py::isinstance<Kmer>(n)) {
        node = dbg.find(n.cast<Kmer>(), true).mappingToFullUnitig();
    } else if (py::isinstance<PyfrostColoredUMap>(n)) {
        node = n.cast<PyfrostColoredUMap>();
    } else {
        throw py::type_error("Unsupported type given, only str, Kmer and UnitigMapping objects supported.");
    }

    if(node.isEmpty || !node.isFullMapping()) {
        throw py::index_error("Node does not exists in the graph");
    }

    return node;
}

template<bool outgoing>
bool EdgeView<outgoing>::isSuccessor(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2) {
    for(auto const& succ : n1.getSuccessors()) {
        if(succ == n2) {
            return true;
        }
    }

    return false;
}

template<bool outgoing>
py::dict EdgeView<outgoing>::get(py::tuple const& edge) {
    py::dict metadata;

    PyfrostColoredUMap n1 = findNode(edge[0]);
    PyfrostColoredUMap n2 = findNode(edge[1]);

    return get(n1, n2);
}

template<bool outgoing>
py::dict EdgeView<outgoing>::get(PyfrostColoredUMap const& n1, PyfrostColoredUMap const& n2) {
    if(!isSuccessor(n1, n2)) {
        throw py::index_error("Target node is not a successor of source node.");
    }

    return makeEdgeDataDict(n1, n2);
}

template<bool outgoing>
EdgeIterator<PyfrostCCDBG::iterator, outgoing> EdgeView<outgoing>::begin() const {
    return {dbg.begin(), dbg.end()};
}

template<bool outgoing>
EdgeIterator<PyfrostCCDBG::iterator, outgoing> EdgeView<outgoing>::end() const {
    return {};
}

template<bool outgoing>
bool EdgeView<outgoing>::contains(const py::tuple& edge) {
    auto n1 = findNode(edge[0]);
    auto n2 = findNode(edge[1]);

    return isSuccessor(n1, n2);
}

template<bool outgoing>
size_t EdgeView<outgoing>::size() const {
    size_t num_edges = 0;
    for(auto const& um : dbg) {
        num_edges += outgoing ? um.getSuccessors().cardinality() : um.getPredecessors().cardinality();
    }

    return num_edges;
}

template<bool outgoing>
class EdgeDataView {
public:
    explicit EdgeDataView(EdgeView<outgoing>& edge_view) : edge_view(edge_view), node_bunch(py::none()),
        data(py::none()), default_value(py::none()) { }

    EdgeDataView(EdgeView<outgoing>& edge_view, py::iterable const& node_bunch) : edge_view(edge_view),
        node_bunch(node_bunch), data(py::none()), default_value(py::none()) { }

    EdgeDataView(EdgeView<outgoing>& edge_view, py::object data, py::object default_value=py::none()) :
        edge_view(edge_view), node_bunch(py::none()), data(std::move(data)), default_value(std::move(default_value)) { }

    EdgeDataView(EdgeView<outgoing>& edge_view, py::iterable const& node_bunch, py::object data,
        py::object default_value=py::none()) :
        edge_view(edge_view), node_bunch(node_bunch), data(std::move(data)),
        default_value(std::move(default_value)) { }

    EdgeDataView(EdgeDataView const& o) = delete;
    EdgeDataView(EdgeDataView&& o) = delete;

    py::object get(py::tuple const& edge) {
        auto n1 = edge_view.findNode(edge[0]);
        auto n2 = edge_view.findNode(edge[1]);

        if(!edge_view.isSuccessor(n1, n2)) {
            throw py::index_error("Not a valid edge");
        }

        auto data_dict = makeEdgeDataDict(n1, n2);

        if(!(data.is_none() || py::isinstance<py::bool_>(data))) {
            auto key = py::cast<std::string>(data);

            if(data_dict.contains(key.c_str())) {
                return data_dict[key.c_str()];
            } else {
                return default_value;
            }
        } else {
            return data_dict;
        }
    }

    bool contains(py::tuple const& edge) {
        if(!edge_view.contains(edge)) {
            return false;
        }

        // Check if data given
        if(edge.size() > 2) {
            auto edge_data = get(edge);
            return edge_data.equal(edge[2]);
        }

        return true;
    }

    bool hasNodeBunch() const {
        return !node_bunch.is_none();
    }

    size_t size() const {
        if(hasNodeBunch()) {
            size_t sum = 0;
            for(auto const& obj : node_bunch) {
                if(!py::isinstance<PyfrostColoredUMap>(obj)) {
                    continue;
                }

                auto um = py::cast<PyfrostColoredUMap>(obj);
                sum += outgoing ? um.getSuccessors().cardinality() : um.getPredecessors().cardinality();
            }

            return sum;
        } else {
            return edge_view.size();
        }
    }

    EdgeWithDataIterator<PyfrostCCDBG::iterator, outgoing> all_begin() const {
        auto edge_iter = EdgeIterator<PyfrostCCDBG::iterator, outgoing>(
            edge_view.dbg.begin(), edge_view.dbg.end());

        return {edge_iter, data, default_value};
    }

    EdgeWithDataIterator<PyfrostCCDBG::iterator, outgoing> all_end() const {
        auto edge_iter = EdgeIterator<PyfrostCCDBG::iterator, outgoing>();
        return {edge_iter, data, default_value};
    }

    EdgeWithDataIterator<NodeBunchIter, outgoing> nbunch_begin() const {
        if(!hasNodeBunch()) {
            throw std::runtime_error("Trying to iterate over node_bunch, but no nodes given!");
        }

        auto edge_iter = EdgeIterator<NodeBunchIter, outgoing>(
            NodeBunchIter(node_bunch.begin()), NodeBunchIter(node_bunch.end()));

        return {edge_iter, data, default_value};
    }

    EdgeWithDataIterator<NodeBunchIter, outgoing> nbunch_end() const {
        auto edge_iter = EdgeIterator<NodeBunchIter, outgoing>();
        return {edge_iter, data, default_value};
    }

private:
    EdgeView<outgoing>& edge_view;
    py::iterable node_bunch;
    py::object data;
    py::object default_value;
};



void define_EdgeView(py::module& m);

}

#endif //PYFROST_EDGEVIEW_H
