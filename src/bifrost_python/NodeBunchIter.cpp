#include "NodeBunchIter.h"
#include "BifrostDiGraph.h"

namespace pyfrost {

NodeBunchIter::NodeBunchIter(py::iterator iter) : wrapped(std::move(iter)) { }

NodeBunchIter::value_type NodeBunchIter::operator*() {
    if(wrapped != py::iterator::sentinel()) {
        return wrapped->cast<NodeBunchIter::value_type>();
    }

    // return empty UnitigMap
    return {};
}

NodeBunchIter::pointer NodeBunchIter::operator->() {
    static const NodeBunchIter::value_type empty;

    if(py::isinstance<NodeBunchIter::value_type>(*wrapped)) {
        auto& um_ref = wrapped->cast<NodeBunchIter::reference>();
        return &um_ref;
    }

    return &empty;
}

NodeBunchIter &NodeBunchIter::operator++() {
    do {
        ++wrapped;
        if(wrapped == py::iterator::sentinel()) {
            break;
        }
    } while(!py::isinstance<NodeBunchIter::value_type>(*wrapped)
            || !wrapped->cast<NodeBunchIter::value_type>().isFullMapping());

    return *this;
}

NodeBunchIter NodeBunchIter::operator++(int) {
    NodeBunchIter tmp(*this);
    operator++();

    return tmp;
}

bool NodeBunchIter::operator==(NodeBunchIter const &o) {
    return wrapped == o.wrapped;
}

bool NodeBunchIter::operator!=(NodeBunchIter const &o) {
    return !operator==(o);
}

}
