#ifndef PYFROST_NEIGHBORS_H
#define PYFROST_NEIGHBORS_H

namespace pyfrost {

template <typename T>
std::vector<T> colorRestrictedSuccessors(T const& unitig, std::unordered_set<size_t> const& allowed_colors) {
    vector<T> neighbors;

    for(auto& succ: unitig.getSuccessors()) {
        auto colorset = succ.getData()->getUnitigColors(succ);

        T first_kmer(succ);
        if(!first_kmer.strand) {
            first_kmer.dist = succ.len - 1;
        }

        first_kmer.len = 1;
        first_kmer.strand = true;

        for(auto const& color: allowed_colors) {
            if(colorset->contains(first_kmer, color)) {
                neighbors.push_back(succ);
                break;
            }
        }
    }

    return neighbors;
}

template <typename T>
std::vector<T> colorRestrictedPredecessors(T const& unitig, std::unordered_set<size_t> const& allowed_colors) {
    vector<T> predecessors;

    for(auto& pred: unitig.getPredecessors()) {
        auto colorset = pred.getData()->getUnitigColors(pred);

        T last_kmer(pred);
        if(last_kmer.strand) {
            last_kmer.dist = last_kmer.len - 1;
        }

        last_kmer.len = 1;
        last_kmer.strand = true;

        for(auto const& color: allowed_colors) {
            if(colorset->contains(last_kmer, color)) {
                predecessors.push_back(pred);
                break;
            }
        }
    }

    return predecessors;
}

}

#endif //PYFROST_NEIGHBORS_H
