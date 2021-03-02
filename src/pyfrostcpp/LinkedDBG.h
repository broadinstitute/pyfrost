#ifndef PYFROST_LINKEDDBG_H
#define PYFROST_LINKEDDBG_H

#include <KmerIterator.hpp>
#include <Kmer.hpp>
#include <robin_hood.h>

#include "pyfrost.h"
#include "Kmer.h"
#include "JunctionTree.h"

namespace pyfrost {

/**
 * The position of a k-mer in Bifrost is always relative to the unitig's forward strand. This function transforms the
 * position to the mapped k-mer's oriented position.
 *
 * @param unitig The UnitigMap object from Bifrost
 * @return oriented position
 */
template<typename U, typename G>
size_t kmer_pos_oriented(UnitigMap<U, G> const& unitig) {
    return unitig.strand
        ? unitig.dist
        : unitig.dist > 0 ? (unitig.size - unitig.getGraph()->getK() - unitig.dist) : 0;
}

using junction_tree_map = robin_hood::unordered_map<Kmer, std::shared_ptr<JunctionTreeNode>>;

template<typename T>
class LinkedDBG {
public:
    explicit LinkedDBG(T* _graph) : graph(_graph) { }
    LinkedDBG(LinkedDBG const& o) = default;
    LinkedDBG(LinkedDBG&& o) noexcept = default;

    LinkedDBG& operator=(LinkedDBG const& o) = default;

    void addLinksFromSequence(std::string const& sequence);
    std::shared_ptr<JunctionTreeNode> getLinks(Kmer const& kmer);

    size_t numTrees() const {
        return junction_trees.size();
    }

    junction_tree_map& getJunctionTrees() {
        return junction_trees;
    }

private:
    std::shared_ptr<JunctionTreeNode> createOrGetTree(Kmer const& kmer);

    T* graph; // non-owning pointer
    junction_tree_map junction_trees;

};

template<typename T>
void LinkedDBG<T>::addLinksFromSequence(const std::string& sequence)
{
    KmerIterator kmer_iter(sequence.c_str()), kmer_end;
    // Keep track through which branch points this sequence is threaded through, to easily update previous junction
    // trees as the threading process progresses.
    auto nodes_to_annotate = std::vector<std::shared_ptr<JunctionTreeNode>>();

    typename T::unitigmap_t prev_unitig;
    for(; kmer_iter != kmer_end; ++kmer_iter) {
        std::pair<Kmer, int> p = *kmer_iter;
        auto umap = graph->find(p.first);

        // Check if k-mer exists in graph
        if(umap.isEmpty) {
            continue;
        }

        auto unitig = umap.mappingToFullUnitig();
        // If any of its successors have an in-degree > 1, then this node should be annotated with a junction tree.
        bool need_annotation = false;
        for(auto& succ : unitig.getSuccessors()) {
            if(succ.getPredecessors().cardinality() > 1) {
                nodes_to_annotate.push_back(createOrGetTree(unitig.getMappedTail()));
                break;
            }
        }

        prev_unitig = unitig;
        Kmer unitig_head = unitig.getMappedHead();

        size_t unitig_len = umap.size - umap.getGraph()->getK() + 1;
        size_t oriented_pos = kmer_pos_oriented(umap);
        size_t diff_to_unitig_end = unitig_len - oriented_pos - 1;

        for(int i = 0; i < diff_to_unitig_end && kmer_iter != kmer_end; ++i) {
            // Taking a shortcut, skip k-mers of sequence that lie within the unitig (i.e. no branches)
            // Don't care if they match or not, to error correct reads
            ++kmer_iter;
        }

        if(kmer_iter == kmer_end) {
            // No branches encountered anymore
            break;
        }

        p = *kmer_iter;
        if(p.first != unitig.getMappedTail()) {
            // Sequence doesn't match with unitig, likely a removed sequencing artifact. There's no unambiguous path
            // that this sequence takes, so we quit
            std::cout << p.first.toString() << " != " << unitig.getMappedTail().toString() << "\n";
            break;
        }

        // We should be at the end of a unitig, so either we have a branch point, or our successor has a in-degree > 1
        size_t next_pos = p.second + Kmer::k;
        if(next_pos >= sequence.length()) {
            // We're already at the end of the sequence, quit
            break;
        }

        if(unitig.getSuccessors().cardinality() > 1) {
            // We have found an unitig with two or more outgoing edges, let's see which path
            // the read takes

            char edge = std::toupper(sequence[next_pos]);
            if(!(edge == 'A' || edge == 'C' || edge == 'G' || edge == 'T')) {
                // Invalid character in sequence, continue, KmerIterator will take care of skipping ahead
                continue;
            }

            Kmer next_kmer = p.first.forwardBase(edge);
            bool found_succ = false;

            // Let's see if there's an successor that matches with the read
            for(auto& succ : umap.getSuccessors()) {
                if(succ.getMappedHead() == next_kmer) {
                    // Add edge choice to each tree
                    std::transform(
                        nodes_to_annotate.begin(), nodes_to_annotate.end(), nodes_to_annotate.begin(),
                        [edge] (std::shared_ptr<JunctionTreeNode> node) {
                            return node->addEdge(edge);
                        });

                    found_succ = true;
                    break;
                }
            }

            if(!found_succ) {
                break;
            }

        }
    }
}

template<typename T>
std::shared_ptr<JunctionTreeNode> LinkedDBG<T>::createOrGetTree(Kmer const& kmer)
{
    auto it = junction_trees.find(kmer);
    if(it == junction_trees.end()) {
        junction_trees.emplace(kmer, std::make_shared<JunctionTreeNode>());
        return junction_trees[kmer];
    } else {
        return it->second;
    }
}

template<typename T>
std::shared_ptr<JunctionTreeNode> LinkedDBG<T>::getLinks(const Kmer& kmer)
{
    return junction_trees.at(kmer);
}

void define_LinkedDBG(py::module& m);

}

#endif //PYFROST_LINKEDDBG_H
