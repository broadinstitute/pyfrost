#include "NodeBunchIter.h"
#include "BifrostDiGraph.h"

namespace pyfrost {

NodeBunchIter::value_type NodeBunchIter::operator*() {
    return wrapped->cast<NodeBunchIter::value_type>();
}

NodeBunchIter &NodeBunchIter::operator++() {
    ++wrapped;

    if(wrapped != py::iterator::sentinel()) {
        auto node = wrapped->cast<NodeBunchIter::reference>();
        while (node.isEmpty || !node.isFullMapping()) {
            ++wrapped;
            if(wrapped == py::iterator::sentinel()) {
                break;
            }

            node = wrapped->cast<NodeBunchIter::reference>();
        }
    }

    return *this;
}

NodeBunchIter NodeBunchIter::operator++(int) {
    NodeBunchIter tmp(*this);
    operator++();

    return tmp;
}

bool NodeBunchIter::operator==(NodeBunchIter const &o) {
    return graph == o.graph && wrapped == o.wrapped;
}

bool NodeBunchIter::operator!=(NodeBunchIter const &o) {
    return !operator==(o);
}

}
