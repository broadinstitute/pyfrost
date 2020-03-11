#include <CompactedDBG.hpp>

#ifndef PYFROST_BIFROSTGRAPH_H
#define PYFROST_BIFROSTGRAPH_H

namespace pyfrost {

class BifrostGraph {
public:
    BifrostGraph() { }
    BifrostGraph(CompactedDBG<>&& dbg) : dbg(std::move(dbg)) { }

private:
    CompactedDBG<> dbg;

};

}


#endif //PYFROST_BIFROSTGRAPH_H
