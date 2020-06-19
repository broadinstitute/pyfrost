#ifndef PYFROST_EDGEVIEW_H
#define PYFROST_EDGEVIEW_H

#include <numeric>
#include <utility>

#include "Pyfrost.h"
#include "NodeIterator.h"

namespace pyfrost {

py::dict makeEdgeDataDict(PyfrostCCDBG& dbg, Kmer const& n1, Kmer const& n2);

// Forward declarations
template<typename T=void, bool outgoing=true> class EdgeWithDataIterator;
template<typename T=void, bool outgoing=true> class EdgeDataView;

/**
 * Iterate over all incoming or outgoing edges for a bunch of nodes. The template parameter T specifies what kind of
 * node container to use, for example the whole graph or a python list.
 *
 * @see NodeIterable
 * @tparam T The node container, if void, uses the whole graph
 * @tparam outgoing Whether to iterate over outgoing edges, else iterates over incoming edges.
 */
template<typename T=void, bool outgoing=true>
class EdgeIterator {
public:
    friend class EdgeWithDataIterator<T, outgoing>;

    using iterator_category = std::input_iterator_tag;
    using value_type = py::tuple;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    using node_iterator_type = typename NodeIterable<T>::iterator_type;
    using edge_iterator_type = PyfrostColoredUMap::neighbor_iterator;

    EdgeIterator(node_iterator_type iter, node_iterator_type end);
    EdgeIterator(EdgeIterator const& o) = default;
    EdgeIterator(EdgeIterator&& o) = default;

    value_type operator*();
    pointer operator->();

    EdgeIterator& operator++();
    EdgeIterator const operator++(int);

    bool operator==(EdgeIterator const& o) const;
    bool operator!=(EdgeIterator const& o) const;

private:
    node_iterator_type node_iter;
    node_iterator_type node_end;

    edge_iterator_type edge_iter;
    edge_iterator_type edge_end;

    value_type current;

    void setCurrentEdge();
};

template<typename T, bool outgoing>
EdgeIterator<T, outgoing>::EdgeIterator(EdgeIterator<T, outgoing>::node_iterator_type iter,
    EdgeIterator<T, outgoing>::node_iterator_type end) : node_iter(iter), node_end(end)
{
    // find a node with outgoing or incoming edges
    while(node_iter != node_end) {
        auto unitig = node_iter.getGraph()->find(*node_iter, true).mappingToFullUnitig();

        if(unitig.isEmpty) {
            ++node_iter;
            continue;
        }

        size_t num = (outgoing ? unitig.getSuccessors().cardinality() :
                      unitig.getPredecessors().cardinality());

        if(num > 0) {
            edge_iter = outgoing ? unitig.getSuccessors().begin() : unitig.getPredecessors().begin();
            edge_end = outgoing ? unitig.getSuccessors().end() : unitig.getPredecessors().end();

            break;
        } else {
            ++node_iter;
        }
    }

    setCurrentEdge();
}

template<typename T, bool outgoing>
void EdgeIterator<T, outgoing>::setCurrentEdge()
{
    if(outgoing) {
        current = py::make_tuple(*node_iter, to_kmer(*edge_iter));
    } else {
        current = py::make_tuple(to_kmer(*edge_iter), *node_iter);
    }
}

template<typename T, bool outgoing>
typename EdgeIterator<T, outgoing>::value_type EdgeIterator<T, outgoing>::operator*()
{
    return current;
}

template<typename T, bool outgoing>
typename EdgeIterator<T, outgoing>::pointer EdgeIterator<T, outgoing>::operator->()
{
    return &current;
}

template<typename T, bool outgoing>
EdgeIterator<T, outgoing>& EdgeIterator<T, outgoing>::operator++()
{
    if(node_iter == node_end) {
        return *this;
    }

    do {
        // Move to next edge if possible
        ++edge_iter;

        if(edge_iter == edge_end) {
            // All edges of current node visited, move to next node
            do {
                ++node_iter;
            } while(is_kmer_empty(*node_iter) && node_iter != node_end);

            if(node_iter == node_end) {
                // Fully done, quit
                break;
            }

            auto unitig = node_iter.getGraph()->find(*node_iter, true).mappingToFullUnitig();
            if(unitig.isEmpty) {
                // Invalid node in the given node iterable. Skip.
                continue;
            }

            if(outgoing) {
                edge_iter = unitig.getSuccessors().begin();
                edge_end = unitig.getSuccessors().end();
            } else {
                edge_iter = unitig.getPredecessors().begin();
                edge_end = unitig.getPredecessors().end();
            }
        }
    } while (edge_iter == edge_end);

    if(node_iter != node_end) {
        setCurrentEdge();
    }

    return *this;
}

template<typename T, bool outgoing>
EdgeIterator<T, outgoing> const EdgeIterator<T, outgoing>::operator++(int)
{
    EdgeIterator<T, outgoing> tmp(*this);
    operator++();

    return tmp;

}

template<typename T, bool outgoing>
bool EdgeIterator<T, outgoing>::operator==(EdgeIterator<T, outgoing> const& o) const {
    if(node_iter == node_end && o.node_iter == o.node_end) {
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
 * @tparam T The node container type, if void uses the whole graph (and thus all nodes/edges).
 * @tparam outgoing
 */
template<typename T, bool outgoing>
class EdgeWithDataIterator {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = py::tuple;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    EdgeWithDataIterator(EdgeIterator<T, outgoing> const& to_wrap, py::object data, py::object default_value) :
        wrapped(to_wrap), data(std::move(data)), default_value(std::move(default_value))
        { }

    EdgeWithDataIterator(EdgeWithDataIterator const& o) = default;
    EdgeWithDataIterator(EdgeWithDataIterator&& o) = default;

    value_type operator*();

    bool operator==(EdgeWithDataIterator const& o) const;
    bool operator!=(EdgeWithDataIterator const& o) const;

    EdgeWithDataIterator& operator++();
    EdgeWithDataIterator const operator++(int);

private:
    EdgeIterator<T, outgoing> wrapped;
    py::object data;
    py::object default_value;

    value_type current;

};

template<typename T, bool outgoing>
typename EdgeWithDataIterator<T, outgoing>::value_type EdgeWithDataIterator<T, outgoing>::operator*() {
    auto edge = *wrapped;

    if(data.is_none() || (py::isinstance<py::bool_>(data) && !data.cast<bool>())) {
        return edge;
    } else {
        auto n1 = py::cast<Kmer>(edge[0]);
        auto n2 = py::cast<Kmer>(edge[1]);
        auto data_dict = makeEdgeDataDict(*(wrapped.node_iter.getGraph()), n1, n2);

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
bool EdgeWithDataIterator<T, outgoing>::operator==(const EdgeWithDataIterator<T, outgoing> &o) const {
    return wrapped == o.wrapped && data.equal(o.data);
}

template<typename T, bool outgoing>
bool EdgeWithDataIterator<T, outgoing>::operator!=(const EdgeWithDataIterator<T, outgoing> &o) const {
    return !operator==(o);
}

template<typename T, bool outgoing>
EdgeWithDataIterator<T, outgoing>& EdgeWithDataIterator<T, outgoing>::operator++() {
    ++wrapped;

    return *this;
}

template<typename T, bool outgoing>
EdgeWithDataIterator<T, outgoing> const EdgeWithDataIterator<T, outgoing>::operator++(int) {
    EdgeWithDataIterator<T, outgoing> tmp(*this);

    operator++();
    return tmp;
}


/**
 * Class to obtain metadata of edges or to iterate over all of them
 * @tparam outgoing Whether to access outgoing edges or the incoming edges
 */
template<bool outgoing>
class EdgeView {
public:
    template<typename U, bool O> friend class EdgeDataView;

    explicit EdgeView(PyfrostCCDBG& dbg);
    EdgeView(EdgeView const& o) = default;
    EdgeView(EdgeView&& o) noexcept = default;

    py::dict get(py::tuple const& edge);
    py::dict get(Kmer const& n1, Kmer const& n2);

    EdgeIterator<void, outgoing> begin() const;
    EdgeIterator<void, outgoing> end() const;

    bool contains(py::tuple const& edge);

    virtual size_t size() const;

    inline PyfrostCCDBG& getGraph() const {
        return dbg;
    }

private:
    Kmer findNode(py::object const& n);
    bool isSuccessor(Kmer const& kmer1, Kmer const& kmer2);

    PyfrostCCDBG& dbg;
};


template<bool outgoing>
EdgeView<outgoing>::EdgeView(PyfrostCCDBG& dbg) : dbg(dbg) { }

template<bool outgoing>
Kmer EdgeView<outgoing>::findNode(py::object const& n) {
    Kmer node = to_kmer(n);

    if(is_kmer_empty(node)) {
        throw py::index_error("Node does not exists in the graph");
    }

    return node;
}

template<bool outgoing>
bool EdgeView<outgoing>::isSuccessor(Kmer const& kmer1, Kmer const& kmer2) {
    auto n1 = dbg.find(kmer1, true);

    for(auto const& succ : n1.getSuccessors()) {
        if(succ.getMappedHead() == kmer2) {
            return true;
        }
    }

    return false;
}

template<bool outgoing>
py::dict EdgeView<outgoing>::get(py::tuple const& edge) {
    py::dict metadata;

    Kmer n1 = findNode(edge[0]);
    Kmer n2 = findNode(edge[1]);

    return get(n1, n2);
}

template<bool outgoing>
py::dict EdgeView<outgoing>::get(Kmer const& n1, Kmer const& n2) {
    if(!isSuccessor(n1, n2)) {
        throw py::index_error("Target node is not a successor of source node.");
    }

    return makeEdgeDataDict(dbg, n1, n2);
}

template<bool outgoing>
EdgeIterator<void, outgoing> EdgeView<outgoing>::begin() const {
    auto iterable = NodeIterable<>(dbg);
    return {iterable.begin(), iterable.end()};
}

template<bool outgoing>
EdgeIterator<void, outgoing> EdgeView<outgoing>::end() const {
    auto iterable = NodeIterable<>(dbg);
    return {iterable.end(), iterable.end()};
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

template<typename T, bool outgoing>
class EdgeDataView {
public:
    EdgeDataView(EdgeView<outgoing>& edge_view, NodeIterable<T> const& iterable) : edge_view(edge_view),
        node_bunch(node_bunch), data(py::none()), default_value(py::none()) { }

    EdgeDataView(EdgeView<outgoing>& edge_view, NodeIterable<T> const& node_bunch, py::object data,
        py::object default_value=py::none()) : edge_view(edge_view), node_bunch(node_bunch), data(std::move(data)),
        default_value(std::move(default_value)) { }

    EdgeDataView(EdgeDataView const& o) = default;
    EdgeDataView(EdgeDataView&& o) noexcept = default;

    py::object get(py::tuple const& edge) {
        if(!edge_view.contains(edge)) {
            throw py::index_error("Not a valid edge");
        }

        auto n1 = edge_view.findNode(edge[0]);
        auto n2 = edge_view.findNode(edge[1]);

        auto data_dict = makeEdgeDataDict(edge_view.dbg, n1, n2);

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

    size_t size() const {
        size_t sum = 0;
        for(auto const& kmer : node_bunch) {
            auto um = node_bunch.getGraph().find(kmer, true);
            sum += outgoing ? um.getSuccessors().cardinality() : um.getPredecessors().cardinality();
        }

        return sum;
    }

    EdgeWithDataIterator<T, outgoing> begin() const {
        auto edge_iter = EdgeIterator<T, outgoing>(node_bunch.begin(), node_bunch.end());
        return {edge_iter, data, default_value};
    }

    EdgeWithDataIterator<T, outgoing> end() const {
        auto edge_iter = EdgeIterator<T, outgoing>(node_bunch.end(), node_bunch.end());
        return {edge_iter, data, default_value};
    }

private:
    EdgeView<outgoing>& edge_view;
    NodeIterable<T> node_bunch;
    py::object data;
    py::object default_value;
};



void define_EdgeView(py::module& m);

template<typename T>
void define_EdgeDataView(py::module& m, std::string const& prefix)
{
    std::string const out_classname = prefix + "OutEdgeDataView";
    py::class_<EdgeDataView<T, true>>(m, out_classname.c_str())
        .def("__getitem__", &EdgeDataView<T, true>::get, py::is_operator())
        .def("__contains__", &EdgeDataView<T, true>::contains, py::is_operator())
        .def("__len__", &EdgeDataView<T, true>::size, py::is_operator())
        .def("__iter__", [] (EdgeDataView<T, true> const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    std::string const in_classname = prefix + "InEdgeDataView";
    py::class_<EdgeDataView<T, false>>(m, in_classname.c_str())
        .def("__getitem__", &EdgeDataView<T, false>::get, py::is_operator())
        .def("__contains__", &EdgeDataView<T, false>::contains, py::is_operator())
        .def("__len__", &EdgeDataView<T, false>::size, py::is_operator())
        .def("__iter__", [] (EdgeDataView<T, false> const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());
}

}

#endif //PYFROST_EDGEVIEW_H
