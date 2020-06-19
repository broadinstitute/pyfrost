#ifndef PYFROST_DEGREEVIEW_H
#define PYFROST_DEGREEVIEW_H

#include <utility>

#include "UnitigDataProxy.h"
#include "NodeIterator.h"


namespace pyfrost {

// Forward declarations
template<typename T=void> class DegreeView;
template<typename T=void> class InDegreeView;
template<typename T=void> class OutDegreeView;

/**
 * Iterate over a bunch of nodes, and return (node, degree) tuples
 *
 * @tparam T The type of the container with nodes (py::list, or PyfrostCCDBG for example).
 */
template<typename T>
class DegreeIterator {
public:
    using node_iterator_type = typename NodeIterable<T>::iterator_type;

    using iterator_category = std::input_iterator_tag;
    using value_type = py::tuple;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    DegreeIterator(DegreeView<T>* view, node_iterator_type wrapped) : view(view), wrapped(wrapped) {
        setCurrentDegree();
    }

    DegreeIterator(DegreeIterator const& o) : view(o.view), wrapped(o.wrapped) {
        setCurrentDegree();
    }

    DegreeIterator(DegreeIterator&& o) noexcept : view(std::move(o.view)), wrapped(std::move(o.wrapped)) {
        setCurrentDegree();
        o.view = nullptr;
    }

    DegreeIterator& operator=(DegreeIterator<T> const& o) {
        if(&o != this) {
            view = o.view;
            wrapped = o.wrapped;
        }

        return *this;
    }

    value_type operator*() {
        if(is_kmer_empty(*wrapped)) {
            return py::make_tuple(py::none(), py::none());
        }

        return py::make_tuple(*wrapped, curr_degree);
    }

    DegreeIterator& operator++() {
        ++wrapped;
        setCurrentDegree();

        return *this;
    }

    DegreeIterator operator++(int) {
        DegreeIterator<T> tmp(*this);
        operator++();

        return tmp;
    }

    bool operator==(DegreeIterator const& o) {
        return view == o.view && wrapped == o.wrapped;
    }

    bool operator!=(DegreeIterator const& o) {
        return !operator==(o);
    }

private:
    DegreeView<T>* view = nullptr;
    node_iterator_type wrapped;
    size_t curr_degree = 0;

    void setCurrentDegree() {
        auto next_node = *wrapped;

        if(!is_kmer_empty(next_node) && view != nullptr) {
            curr_degree = view->getDegree(next_node);
        }
    }
};

/**
 * Obtain the degree for a bunch of nodes
 *
 * The type of container with nodes can specified using the template param. By default iterates over all nodes in the
 * graph.
 * @tparam T
 */
template<typename T>
class DegreeView {
public:
    using iterator_type = DegreeIterator<T>;
    friend class DegreeIterator<T>;

    explicit DegreeView(NodeIterable<T> const& node_bunch) : node_bunch(node_bunch) { }

    virtual ~DegreeView() = default;

    iterator_type begin() {
        return iterator_type(this, node_bunch.begin());
    }

    iterator_type end() {
        return iterator_type(this, node_bunch.end());
    }

    virtual size_t getDegree(Kmer const& kmer) const {
        size_t degree = 0;
        auto unitig = node_bunch.getGraph().find(kmer, true);

        if(unitig.isEmpty) {
            throw py::key_error("Node not found.");
        }

        for(auto const& p : unitig.getPredecessors()) {
            ++degree;
        }

        for(auto const& s : unitig.getSuccessors()) {
            ++degree;
        }

        return degree;
    }

    size_t getDegree(char const* kmer) const {
        return getDegree(Kmer(kmer));
    }

    bool contains(Kmer const& kmer) const {
        auto unitig = node_bunch.getGraph().find(kmer, true);
        return !unitig.isEmpty;
    }

    PyfrostCCDBG& getGraph() {
        return node_bunch.getGraph();
    }

protected:
    NodeIterable<T> node_bunch;
};


template<typename T>
class InDegreeView : public DegreeView<T> {
public:
    explicit InDegreeView(NodeIterable<T> const& node_bunch) : DegreeView<T>(node_bunch) { }
    ~InDegreeView() override = default;

    size_t getDegree(Kmer const& kmer) const override {
        size_t degree = 0;
        auto unitig = this->node_bunch.getGraph().find(kmer, true);

        if(unitig.isEmpty) {
            throw py::key_error("Node not found.");
        }

        for(auto const& p : unitig.getPredecessors()) {
            ++degree;
        }

        return degree;
    }
};

template<typename T>
class OutDegreeView : public DegreeView<T> {
public:
    explicit OutDegreeView(NodeIterable<T> const& node_bunch) : DegreeView<T>(node_bunch) { }
    ~OutDegreeView() override = default;

    size_t getDegree(Kmer const& kmer) const override {
        size_t degree = 0;
        auto unitig = this->node_bunch.getGraph().find(kmer, true);

        if(unitig.isEmpty) {
            throw py::key_error("Node not found.");
        }

        for(auto const& s : unitig.getSuccessors()) {
            ++degree;
        }

        return degree;
    }
};

template<typename T>
void define_DegreeView(py::module& m, std::string const& class_prefix) {

    std::string const degreeview_name = class_prefix + "DegreeView";
    auto py_DegreeView = py::class_<DegreeView<T>>(m, degreeview_name.c_str())
        .def("__getitem__", [] (DegreeView<T>& self, Kmer const& kmer) {
            return self.getDegree(kmer);
        })
        .def("__getitem__", [] (DegreeView<T>& self, char const* kmer) {
            return self.getDegree(kmer);
        })

        .def("__call__", [] (DegreeView<T>& self) {
            return self;
        })

        // Call with just single node as argument
        .def("__call__", [] (DegreeView<T>& self, Kmer const& kmer, py::object const&) {
            // Weight is ignored
            return self.getDegree(kmer);
        }, py::arg("nbunch"), py::arg("weight") = py::none())
        .def("__call__", [] (DegreeView<T>& self, char const* kmer, py::object const&) {
            return self.getDegree(kmer);
        }, py::arg("nbunch"), py::arg("weight") = py::none())

        // Call with a bunch of nodes as arguments. Return a new degree view
        .def("__call__", [] (DegreeView<T>& self, py::iterable const& nbunch, py::object const&) {
            return DegreeView<py::iterable>(NodeIterable<py::iterable>(self.getGraph(), nbunch));
        }, py::arg("nbunch"), py::arg("weight") = py::none())

        // Iterate over nodes, either in the container or over the whole graph.
        .def("__iter__", [] (DegreeView<T>& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>());

    auto Iterable = py::module::import("collections.abc").attr("Iterable");
    py_DegreeView.attr("__bases__") = py::make_tuple(Iterable).attr("__add__")(py_DegreeView.attr("__bases__"));

    std::string const indegreeview_name = class_prefix + "InDegreeView";
    py::class_<InDegreeView<T>, DegreeView<T>>(m, indegreeview_name.c_str())
        .def("__getitem__", [] (InDegreeView<T>& self, Kmer const& kmer) {
            return self.getDegree(kmer);
        })
        .def("__getitem__", [] (DegreeView<T>& self, char const* kmer) {
            return self.getDegree(kmer);
        })

        // Call with a bunch of nodes as arguments. Make sure to return a new indegree view
        .def("__call__", [] (DegreeView<T>& self, py::iterable const& nbunch, py::object const&) {
             return InDegreeView<py::iterable>(NodeIterable<py::iterable>(self.getGraph(), nbunch));
         }, py::arg("nbunch"), py::arg("weight") = py::none());

    std::string const outdegreeview_name = class_prefix + "OutDegreeView";
    py::class_<OutDegreeView<T>, DegreeView<T>>(m, outdegreeview_name.c_str())
        .def("__getitem__", [] (OutDegreeView<T>& self, Kmer const& kmer) {
            return self.getDegree(kmer);
        })
        .def("__getitem__", [] (DegreeView<T>& self, char const* kmer) {
            return self.getDegree(kmer);
        })

        // Call with a bunch of nodes as arguments. Make sure to return a new indegree view
        .def("__call__", [] (DegreeView<T>& self, py::iterable const& nbunch, py::object const&) {
             return OutDegreeView<py::iterable>(NodeIterable<py::iterable>(self.getGraph(), nbunch));
         }, py::arg("nbunch"), py::arg("weight") = py::none());
}

}

#endif //PYFROST_DEGREEVIEW_H
