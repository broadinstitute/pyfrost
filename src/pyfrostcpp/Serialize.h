#ifndef PYFROST_SERIALIZE_H
#define PYFROST_SERIALIZE_H

#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <robin_hood.h>
#include <tl/optional.hpp>

// Cereal serialization support for robin_hood::unordered_map
namespace robin_hood {

template <typename Archive, typename K, typename V, typename C, typename A>
void save(Archive& ar, unordered_map<K, V, C, A> const& map) {
    ar(map.size());
    for(auto const& e: map) {
        ar(e.first, e.second);
    }
}

template <typename Archive, typename K, typename V, typename C, typename A>
void load(Archive& ar, unordered_map<K, V, C, A>& map) {
    map.clear();

    size_t num_entries;
    ar(num_entries);
    map.reserve(num_entries);

    for(int i = 0; i < num_entries; ++i) {
        K key;
        V value;
        ar(key);
        ar(value);

        map.emplace(std::move(key), std::move(value));
    }
}

}

namespace tl {

template<typename Archive, typename T>
void save(Archive& ar, optional<T> const& value) {
    if(value) {
        ar(true, value.value());
    } else {
        ar(false);
    }
}

template<typename Archive, typename T>
void load(Archive& ar, optional<T>& opt_value) {
    bool has_value;
    ar(has_value);
    if(has_value) {
        T value;
        ar(value);
        opt_value = value;
    } else {
        opt_value = tl::nullopt;
    }
}

}

#endif //PYFROST_SERIALIZE_H
