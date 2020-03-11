#include <ColoredCDBG.hpp>

#ifndef PYFROST_COLOREDCDBGRAPH_H
#define PYFROST_COLOREDCDBGRAPH_H

namespace pyfrost {

template<typename U=void>
class ColoredBifrostGraph {
public:
    ColoredBifrostGraph() { }

    ColoredBifrostGraph(ColoredCDBG<U>&& dbg) {
        ccdbg = std::move(dbg);
    }

    int num_kmers() const {
        return ccdbg.nbKmers();
    }

private:
    ColoredCDBG<U> ccdbg;
};

}


#endif //PYFROST_COLOREDCDBGRAPH_H
