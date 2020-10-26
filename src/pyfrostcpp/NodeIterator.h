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
 *
 * Can optionally also include reverse complement of nodes/kmers
 */
template<typename T>
class NodeIterator {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = Kmer;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    NodeIterator(PyfrostCCDBG* dbg, T iter, bool _with_rev_compl=false) :
        dbg(dbg), wrapped(iter), is_rev_compl(false), with_rev_compl(_with_rev_compl)
    {
        setCurrentKmer();
    }

    NodeIterator(NodeIterator const& o) : dbg(o.dbg), wrapped(o.wrapped), is_rev_compl(o.is_rev_compl),
        with_rev_compl(o.with_rev_compl) {
        setCurrentKmer();
    };

    NodeIterator(NodeIterator&& o) noexcept : dbg(o.dbg), wrapped(std::move(o.wrapped)), is_rev_compl(o.is_rev_compl),
        with_rev_compl(o.with_rev_compl)
    {
        setCurrentKmer();
        o.dbg = nullptr;
    }

    NodeIterator& operator=(NodeIterator<T> const& o) {
        if(&o != this) {
            dbg = o.dbg;
            wrapped = o.wrapped;
            is_rev_compl = o.is_rev_compl;
            with_rev_compl = o.with_rev_compl;

            setCurrentKmer();
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
        if(with_rev_compl) {
            if(is_rev_compl) {
                ++wrapped;
                is_rev_compl = false;
            } else {
                is_rev_compl = true;
            }
        } else {
            ++wrapped;
        }

        setCurrentKmer();

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
    bool is_rev_compl;
    bool check_valid;
    bool with_rev_compl;

    void setCurrentKmer() {
        current = to_kmer(*wrapped, with_rev_compl && is_rev_compl);
    }

};


}



#endif //PYFROST_NODEITERATOR_H
