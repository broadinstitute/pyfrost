#include <CompactedDBG.hpp>

#ifndef PYFROST_NODEVIEW_H
#define PYFROST_NODEVIEW_H

namespace pyfrost {

template <typename U=void>
class NodeView {
public:
    NodeView() : dbg(nullptr) { }
    NodeView(CompactedDBG<>& dbg) : dbg(&dbg) { }



private:
    CompactedDBG<>* dbg;
};

}

#endif //PYFROST_NODEVIEW_H
