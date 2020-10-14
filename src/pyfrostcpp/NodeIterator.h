#ifndef PYFROST_NODEITERATOR_H
#define PYFROST_NODEITERATOR_H

#include "pyfrost.h"
#include "Kmer.h"
#include "NodeDataDict.h"

std::ostream& operator<<(std::ostream& o, pyfrost::PyfrostColoredUMap const& u);

namespace pyfrost {

/**
 * Iterate over a collection of nodes (unitigs), where each node is represented by its head Kmer.
 *
 * This class is templated, so you can use any kind of iterator as source for nodes. Each item of the wrapped
 * iterator is converted to a Kmer object using `to_kmer`.
 */
template<typename T>
class NodeIterator {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = Kmer;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    NodeIterator(PyfrostCCDBG* dbg, T iter, bool _check_valid=true) :
        dbg(dbg), wrapped(iter), check_valid(_check_valid)
    {
        setCurrentKmer();
    }

    NodeIterator(NodeIterator const& o) : dbg(o.dbg), wrapped(o.wrapped), check_valid(o.check_valid) {
        setCurrentKmer();
    };

    NodeIterator(NodeIterator&& o) noexcept : dbg(o.dbg), wrapped(std::move(o.wrapped)), check_valid(o.check_valid) {
        setCurrentKmer();
        o.dbg = nullptr;
    }

    NodeIterator& operator=(NodeIterator<T> const& o) {
        if(&o != this) {
            dbg = o.dbg;
            wrapped = o.wrapped;
        }

        return *this;
    }

    value_type operator*() {
        return current;
    }

    pointer operator->() {
        return &current;
    }

    /**
     * Move wrapped iterator to the next item, and attempt to convert it to a Kmer object. If resulting k-mer is not
     * a unitig, return an empty k-mer object.
     */
    NodeIterator& operator++() {
        if(check_valid) {
            // Move iterator along. If we encounter an invalid value, continue moving, until the iterator doesn't change
            // anymore.
            T prev;
            do {
                prev = wrapped;
                ++wrapped;

                setCurrentKmer();
            } while (is_kmer_empty(current) && prev != wrapped);

        } else {
            ++wrapped;
            setCurrentKmer();
        }

        return *this;
    }

    NodeIterator const operator++(int) {
        NodeIterator<T> const tmp(*this);
        operator++();

        return tmp;
    }

    bool operator==(NodeIterator const& o) const {
        return wrapped == o.wrapped;
    }

    bool operator!=(NodeIterator const& o) const {
        return !operator==(o);
    }

    PyfrostCCDBG* getGraph() {
        if(dbg == nullptr) {
            throw std::runtime_error("Trying to obtain non-existent graph from NodeIterator.");
        }

        return dbg;
    }

private:
    PyfrostCCDBG* dbg = nullptr;
    T wrapped;
    Kmer current;
    bool check_valid;

    void setCurrentKmer() {
        current = to_kmer(*wrapped);

        if(check_valid) {
            if(dbg != nullptr) {
                auto unitig = dbg->find(current, true);

                if(unitig.isEmpty) {
                    current.set_empty();
                }
            } else {
                current.set_empty();
            }
        }
    }

};

/**
 * This class provides a wrapper around any container with nodes, for easy access to `NodeIterator` objects.
 *
 * Example containers are any python iterables, or some C++ container like vector.
 *
 * @tparam T container typename with nodes
 */
template<typename T=void>
class NodeIterable {
private:
    PyfrostCCDBG& dbg;
    T iterable;

public:
    using iterator_inner_type = decltype(iterable.begin());
    using iterator_type = NodeIterator<iterator_inner_type>;

    NodeIterable(PyfrostCCDBG& dbg, T iterable) : dbg(dbg), iterable(iterable) { }

    NodeIterable(NodeIterable<T> const& o) = default;
    NodeIterable(NodeIterable<T>&& o) noexcept = default;

    iterator_type begin() const {
        return iterator_type(&dbg, iterable.begin());
    }

    iterator_type end() const {
        return iterator_type(&dbg, iterable.end());
    }

    PyfrostCCDBG& getGraph() const {
        return dbg;
    }
};

/**
 * Specialization for NodeIterator without any specific container. Will iterate over all nodes in the graph.
 */
template<>
class NodeIterable<void> {
private:
    PyfrostCCDBG& dbg;

public:
    using iterator_inner_type = decltype(dbg.begin());
    using iterator_type = NodeIterator<iterator_inner_type>;

    explicit NodeIterable(PyfrostCCDBG& dbg) : dbg(dbg) { }

    NodeIterable(NodeIterable<void> const& o) = default;
    NodeIterable(NodeIterable<void>&& o) noexcept = default;

    iterator_type begin() const {
        return {&dbg, dbg.begin(), false};
    }

    iterator_type end() const {
        return {&dbg, dbg.end(), false};
    }

    PyfrostCCDBG& getGraph() const {
        return dbg;
    }
};


}



#endif //PYFROST_NODEITERATOR_H
