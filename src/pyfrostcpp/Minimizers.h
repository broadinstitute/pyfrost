#ifndef PYFROST_MINIMIZERS_H
#define PYFROST_MINIMIZERS_H

#include "pyfrost.h"

#include <minHashIterator.hpp>
#include <pybind11/pybind11.h>

namespace pyfrost {

class MinHashIterWrapper {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = std::vector<minHashResult>;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    MinHashIterWrapper() { }
    MinHashIterWrapper(MinHashIterWrapper const& o) = default;
    MinHashIterWrapper(MinHashIterWrapper&& o) = default;
    MinHashIterWrapper(minHashIterator<RepHash> const& to_wrap) : wrapped_iter(to_wrap) {
        fillResults();
    }

    MinHashIterWrapper& operator++() {
        ++wrapped_iter;
        fillResults();

        return *this;
    }

    MinHashIterWrapper operator++(int) {
        MinHashIterWrapper tmp(*this);

        operator++();
        return tmp;
    }

    reference operator*() {
        return results;
    }

    pointer operator->() {
        return &results;
    }

    bool operator==(MinHashIterWrapper const& o) {
        return wrapped_iter == o.wrapped_iter;
    }

    bool operator!=(MinHashIterWrapper const& o) {
        return !operator==(o);
    }

private:
    minHashIterator<RepHash> wrapped_iter;
    std::vector<minHashResult> results;

    void fillResults() {
        results.clear();

        if(wrapped_iter.invalid) {
            return;
        }

        minHashResultIterator<RepHash> result_iter = *wrapped_iter, result_iter_end;
        for(; result_iter != result_iter_end; ++result_iter) {
            results.push_back(*result_iter);
        }
    }
};

class MinHashIterator {
public:
    explicit MinHashIterator(std::string const& _str) : str(_str), k(Kmer::k), g(Minimizer::g) { }
    MinHashIterator(std::string const& _str, size_t _k) : str(_str), k(_k) {
        setKG(_k);
    }

    MinHashIterator(std::string const& _str, size_t _k, size_t _g) : str(_str), k(_k), g(_g) {
        setKG(_k, _g);
    }

    MinHashIterWrapper begin() const {
        return {minHashIterator<RepHash>(str.c_str(), str.length(), k, g, RepHash(), false)};
    }

    MinHashIterWrapper end() const {
        return {};
    }

private:
    std::string str;
    size_t k;
    size_t g;

};


void define_Minimizer(py::module& m);
void define_MinHashIterator(py::module& m);
void define_MinHashResult(py::module& m);

}

#endif
