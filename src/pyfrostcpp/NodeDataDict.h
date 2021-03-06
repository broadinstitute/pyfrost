#ifndef PYFROST_NODEDATADICT_H
#define PYFROST_NODEDATADICT_H

#include <unordered_map>

#include "pyfrost.h"

namespace pyfrost {

void define_NodeDataDict(py::module&);

enum class UnitigMetaKeys : uint8_t {
    LENGTH,
    POS,
    STRAND,
    HEAD,
    TAIL,
    MAPPED_SEQUENCE,
    UNITIG_SEQUENCE,
    UNITIG_LENGTH,
    IS_FULL_MAPPING,
    COLORS,
    NONE
};

class UnitigDataKeyIterator;

class NodeDataDict {
public:
    friend class UnitigDataKeyIterator;

    explicit NodeDataDict(PyfrostColoredUMap const& _unitig) : unitig(_unitig) {
        if(unitig.isEmpty) {
            throw std::runtime_error("Trying to construct NodeDataDict for non-existent unitig.");
        }
    }
    NodeDataDict(NodeDataDict const& o) = default;
    NodeDataDict(NodeDataDict&& o) = default;

    NodeDataDict& operator=(NodeDataDict const& o) = default;

    bool contains(std::string const& key) const;
    py::object getData(std::string const& key) const;
    void setData(std::string const& key, py::object const& value);
    void delData(py::str const& key);

    size_t size() const;

    std::string mappedSequence() const;
    size_t mappedSequenceLength() const;

    std::string unitigSequence() const;
    size_t unitigLength() const;

    Kmer unitigHead() const;
    Kmer unitigTail() const;
    Kmer unitigRepresentative() const;

    inline PyfrostColoredUMap mappingToFullUnitig() const {
        return unitig.mappingToFullUnitig();
    }

    UnitigDataKeyIterator begin() const;
    UnitigDataKeyIterator end() const;

private:
    py::dict& getDataDict() const {
        auto data_ptr = unitig.getData()->getData(unitig);

        if(data_ptr == nullptr) {
            stringstream error_msg;
            error_msg << "Could not obtain unitig data? nullptr. Head: ";
            error_msg << unitig.getUnitigHead().toString();
            error_msg << " isEmpty: " << unitig.isEmpty;
            throw std::runtime_error(error_msg.str());
        }

        return data_ptr->getDict();
    }

    UnitigMetaKeys getMetaKey(std::string const& key) const;

    PyfrostColoredUMap unitig;
    static const std::unordered_map<std::string, UnitigMetaKeys> hardcoded_keys;
};

/**
 * This iterator wraps both the hardcoded keys and the keys in the user dictionary.
 */
class UnitigDataKeyIterator {
private:
    using map_iter = std::unordered_map<std::string, UnitigMetaKeys>::const_iterator;
    map_iter it1;
    map_iter it1_end;
    py::detail::dict_iterator it2;

public:
    using iterator_category = std::input_iterator_tag;
    using value_type = py::str;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    UnitigDataKeyIterator(map_iter it1, map_iter it1_end, py::detail::dict_iterator it2) :
        it1(it1), it1_end(it1_end), it2(it2) { }

    UnitigDataKeyIterator(UnitigDataKeyIterator const& o) = default;
    UnitigDataKeyIterator(UnitigDataKeyIterator&& o) = default;

    value_type operator*() {
        if(it1 == it1_end) {
            // All hardcoded keys done, now return the keys in the user dict
            return it2->first.cast<value_type>();
        } else {
            return it1->first;
        }
    }

    UnitigDataKeyIterator& operator++() {
        if(it1 == it1_end) {
            ++it2;
        } else {
            ++it1;
        }

        return *this;
    }

    UnitigDataKeyIterator operator++(int) {
        const UnitigDataKeyIterator tmp(*this);

        operator++();

        return tmp;
    }

    bool operator==(UnitigDataKeyIterator const& o) const {
        if(&o == this) {
            return true;
        }

        return it1 == o.it1 && it1_end == o.it1_end && it2 == o.it2;
    }

    bool operator!=(UnitigDataKeyIterator const& o) {
        return !operator==(o);
    }

};

}

#endif //PYFROST_NODEDATADICT_H
