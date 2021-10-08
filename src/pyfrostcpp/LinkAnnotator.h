#ifndef PYFROST_LINKANNOTATOR_H
#define PYFROST_LINKANNOTATOR_H

#include "pyfrost.h"
#include "JunctionTree.h"
#include "LinkDB.h"

#include <condition_variable>

namespace pyfrost {

struct MappingResult {
    MappingResult() : mapping_start(0), mapping_end(0) { }
    MappingResult(MappingResult const& o) = default;
    MappingResult(MappingResult&& o) = default;

    /// Path through the graph the sequence took
    vector<Kmer> path;

    /// First position in the sequence that has a k-mer that maps to the graph
    size_t mapping_start;

    /// Last position in the sequence that has a k-mer that maps to the graph
    size_t mapping_end;

    /// Junction choices from this link
    string junctions;

    /// For each k-mer in the sequence, whether it was found in the graph or not
    std::vector<uint8_t> matches;

    /// For each unitig this sequence traversed, count how often we encountered it
    std::unordered_map<Kmer, size_t> unitig_visits;

    Kmer start_unitig() const {
        if(path.empty()) {
            Kmer empty;
            empty.set_empty();
            return empty;
        } else {
            return path[0];
        }
    }

    Kmer end_unitig() const {
        if(path.empty()) {
            Kmer empty;
            empty.set_empty();
            return empty;
        } else {
            return path[path.size()-1];
        }
    }
};


/**
 * The position of a k-mer in Bifrost is always relative to the unitig's forward strand. This function transforms the
 * position to the mapped k-mer's oriented position.
 *
 * @param unitig The UnitigMap object from Bifrost
 * @return oriented position
 */
template<typename U, typename G>
size_t kmerPosOriented(UnitigMap<U, G> const& unitig) {
    return unitig.strand
           ? unitig.dist
           : unitig.size - unitig.getGraph()->getK() - unitig.dist;
}


template<typename T>
class LinkAnnotator {
public:
    LinkAnnotator(T* _graph, LinkDB* _db) : graph(_graph), db(_db) { }
    LinkAnnotator(LinkAnnotator const& o) = default;
    LinkAnnotator(LinkAnnotator&& o) noexcept = default;
    virtual ~LinkAnnotator() = default;

    /**
     * Thread the given sequence through the graph and annotate specific nodes with the junctions taken.
     *
     * @param seq the sequence to thread through the graph
     */
    virtual MappingResult addLinksFromSequence(std::string const& seq);

    /**
     * Thread the given sequence through the graph and annotate specific nodes with the junctions taken.
     *
     * If calling `addLinksFromSequence` multiple times, by default each sequence will be considered
     * independently, and previously annotated nodes will not be considered for the new sequence. In case of paired-end
     * reads, however, when processing the second read in a pair it might be worthwhile to continue adding links to
     * nodes annotated by the first read. In that case you can set `keep_nodes` to true.
     *
     * @param seq the sequence to thread through the graph
     * @param keep_nodes Whether to keep the nodes to annotate from a previous call
     */
    virtual MappingResult addLinksFromSequence(std::string const& seq, bool keep_nodes);

    /**
     * Follow the given path through the graph, and add the junction choices to the current junction trees.
     *
     * This is useful when mapping paired reads in two passes. In the first pass, you generate links from reads
     * ignoring read-pair information. Then in a second pass, you separately generate links again, but now utilizing
     * the links generated in the first pass to find a link supported path from the end unitig of read 1, to the start
     * unitig of read 2. In that way you can continue adding links to the junction trees generated from read 1, with
     * junction choices from read 2.
     */
    virtual void addLinksFromPath(vector<Kmer> const& path);

protected:
    /**
     * This function assesses whether a given node (unitig) needs link annotations.
     *
     * By default, any node which has at least one successor with in degree > 1 or a predecessor with out-degree > 1
     * qualifies as a node for link annotations.
     *
     * See paper:
     * Turner, Isaac, et al. "Integrating long-range connectivity information into de Bruijn graphs."
     * Bioinformatics 34.15 (2018): 2556-2565.
     *
     * @tparam T Type of graph: CompactedDBG, ColoredCDBG, PyfrostCCDBG, ...
     * @param unitig the node to check
     * @return True if nodes needs annotations
     */
    virtual bool nodeNeedsAnnotation(typename T::unitigmap_t const& unitig);

    virtual void resetNodesToAnnotate() {
        nodes_to_annotate.clear();
    }

    /**
     * Locate a k-mer in the graph. Can be overriden to ensure only k-mers with a color are valid.
     */
    virtual typename T::unitigmap_t findKmer(Kmer kmer) {
        return graph->find(kmer);
    }

    /**
     * Returns number of successors of a given node, can be overridden for e.g. color specific successors.
     */
    virtual size_t numSucessors(typename T::unitigmap_t const& unitig) {
        return unitig.getSuccessors().cardinality();
    }

    /**
     * Returns number of predecessors of a given node, can be overridden for e.g. color specific predecessors.
     */
    virtual size_t numPredecessors(typename T::unitigmap_t const& unitig) {
        return unitig.getPredecessors().cardinality();
    }

    T* graph;
    LinkDB* db;
    std::vector<JunctionTreeNode*> nodes_to_annotate;
};

template<typename T>
bool LinkAnnotator<T>::nodeNeedsAnnotation(typename T::unitigmap_t const& unitig) {
    for(auto& succ : unitig.getSuccessors()) {
        if(numPredecessors(succ) > 1) {
            return true;
        }
    }

    for(auto& pred : unitig.getPredecessors()) {
        if(numSucessors(pred) > 1) {
            return true;
        }
    }

    return false;
}

template<typename T>
MappingResult LinkAnnotator<T>::addLinksFromSequence(string const& seq) {
    return LinkAnnotator<T>::addLinksFromSequence(seq, false);
}

template<typename T>
MappingResult LinkAnnotator<T>::addLinksFromSequence(string const& seq, bool keep_nodes) {
    if(graph == nullptr || db == nullptr) {
        throw std::runtime_error("Graph or LinkDB pointer expired!");
    }

    if(!keep_nodes) {
        resetNodesToAnnotate();
    }

    bool first_unitig_found = false;
    MappingResult mapping;
    mapping.matches.resize(seq.length() - Kmer::k + 1, 0);

    KmerIterator kmer_iter(seq.c_str()), kmer_end;
    for(; kmer_iter != kmer_end; ++kmer_iter) {
        Kmer kmer;
        size_t pos;
        tie(kmer, pos) = *kmer_iter;

        auto umap = findKmer(kmer);
        if(umap.isEmpty) {
            if(first_unitig_found){
                break;
            } else {
                // Allow for "clipping" of the sequence if we haven't found a unitig yet.
                continue;
            }
        }

        auto unitig = umap.mappingToFullUnitig();
        auto unitig_kmer = unitig.getMappedHead();

        if(!first_unitig_found) {
            // First k-mer that is present in the graph
            mapping.mapping_start = pos;
            first_unitig_found = true;
        }

        mapping.path.push_back(unitig_kmer);
        mapping.matches[pos] = true;
        mapping.mapping_end = pos;

        auto it = mapping.unitig_visits.find(unitig_kmer);
        if(it != mapping.unitig_visits.end()) {
            ++(it->second);
        } else {
            mapping.unitig_visits[unitig_kmer] = 1;
        }

        if(nodeNeedsAnnotation(unitig)) {
            nodes_to_annotate.push_back(&db->createOrGetTree(unitig.getMappedTail()));
        }

        // Move to the end of the unitig, and by definition we will not encounter any branch points.
        size_t unitig_len = unitig.size - graph->getK() + 1; // unitig length in num k-mers
        size_t oriented_pos = kmerPosOriented(umap);
        size_t diff_to_unitig_end = unitig_len - oriented_pos - 1;

        for(int i = 0; i < diff_to_unitig_end && kmer_iter != kmer_end; ++i) {
            // To error correct reads, we don't care if the kmers of the unitig don't match the k-mers of the given
            // sequence, as long as it ends up on the same unitig again (see code below for loop). We do keep track
            // of mismatches for analysis purposes.
            ++kmer_iter;

            if(kmer_iter != kmer_end) {
                tie(kmer, pos) = *kmer_iter;

                size_t unitig_pos = umap.dist + (umap.strand ? i + 1 : -i - 1);
                unitig_kmer = unitig.getMappedKmer(unitig_pos);
                bool match = kmer == unitig_kmer;
                mapping.matches[pos] = match;
                if(match) {
                    mapping.mapping_end = pos;
                }
            }
        }

        if(kmer_iter == kmer_end) {
            // End of sequence, but still on the same unitig, so no branches encountered
            break;
        }

        if(kmer != unitig.getMappedTail()) {
            // The k-mer doesn't match the unitig tail anymore. The given sequence (a read) likely contained an error
            // that was removed from the graph, and now there's no unambiguous path anymore, so we quit.
            break;
        }

        size_t edge_pos = pos + Kmer::k;
        if(edge_pos >= seq.length()) {
            // At the end of the sequence, so no edge to inspect
            break;
        }

        if(numSucessors(unitig) > 1) {
            char edge = toupper(seq[edge_pos]);
            if(!(edge == 'A' || edge == 'C' || edge == 'G' || edge == 'T')) {
                // Invalid sequence, quit
                break;
            }

            Kmer succ_kmer = kmer.forwardBase(edge);
            bool found_succ = false;

            // Let's see if there's a successor that matches with the sequence
            for(auto& succ : unitig.getSuccessors()) {
                if(succ.getMappedHead() == succ_kmer) {
                    // Add edge choice to each tree
                    std::transform(
                        nodes_to_annotate.begin(), nodes_to_annotate.end(), nodes_to_annotate.begin(),
                        [edge] (JunctionTreeNode* node) -> JunctionTreeNode* {
                            return &node->addEdge(edge);
                        });

                    found_succ = true;
                    break;
                }
            }

            // Couldn't find a valid successor in the graph, abort
            if(!found_succ) {
                break;
            }
        }
    }

    if(!nodes_to_annotate.empty()) {
        mapping.junctions = nodes_to_annotate[0]->getJunctionChoices();
    }

    return mapping;
}

template<typename T>
void LinkAnnotator<T>::addLinksFromPath(vector<Kmer> const& path) {
    if(graph == nullptr || db == nullptr) {
        throw std::runtime_error("Graph or LinkDB pointer expired!");
    }

    size_t i = 0;
    for(auto const& kmer : path) {
        auto unitig = graph->find(kmer, true).mappingToFullUnitig();
        if(unitig.isEmpty) {
            throw std::runtime_error("Invalid path! Kmer " + kmer.toString() + " is not a unitig.");
        }

        // Don't mark the start and end of the path as nodes to annotate, because those nodes will be already be marked
        // as to be annotated by `addLinksFromSequence`.
        if(i > 0 && i < path.size() - 1) {
            if(nodeNeedsAnnotation(unitig)) {
                nodes_to_annotate.push_back(&db->createOrGetTree(unitig.getMappedTail()));
            }
        }

        // Follow path, and see which edge is taken at junctions
        if(i < path.size() - 1 && numSucessors(unitig) > 1) {
            char edge = path[i+1].getChar(Kmer::k-1);
            std::transform(
                nodes_to_annotate.begin(), nodes_to_annotate.end(), nodes_to_annotate.begin(),
                [edge] (JunctionTreeNode* node) -> JunctionTreeNode* {
                    return &node->addEdge(edge);
                });
        }

        ++i;
    }
}

/**
 * Annotator class that only follows a given color. I.e., only junctions for a given color are annotated.
 * @tparam T
 */
template<typename T>
class ColorAssociatedAnnotator : public LinkAnnotator<T> {
public:
    ColorAssociatedAnnotator(T* _graph, LinkDB* _linkdb, size_t _color)
        : LinkAnnotator<T>(_graph, _linkdb), color(_color)
    { }
    ColorAssociatedAnnotator(ColorAssociatedAnnotator const& o) = default;
    ColorAssociatedAnnotator(ColorAssociatedAnnotator&& o) noexcept = default;

    virtual ~ColorAssociatedAnnotator() = default;

protected:
    typename T::unitigmap_t findKmer(Kmer kmer) override {
        auto umap = this->graph->find(kmer);
        if(umap.isEmpty) {
            return umap;
        }

        auto colorset = umap.getData()->getUnitigColors(umap);

        if(colorset->contains(umap, color)) {
            return umap;
        } else {
            typename T::unitigmap_t empty;
            return empty;
        }
    }

    size_t numSucessors(typename T::unitigmap_t const& unitig) override {
        size_t num_with_color = 0;
        for(auto& succ : unitig.getSuccessors()) {
            auto colorset = succ.getData()->getUnitigColors(succ);

            // Make mapping only include first or last k-mer (depending on whether we are on the reverse complement
            // or not)
            typename T::unitigmap_t first_kmer(succ);
            if(!first_kmer.strand) {
                first_kmer.dist = succ.len - 1;
            }

            first_kmer.len = 1;
            first_kmer.strand = true;

            if(colorset->contains(first_kmer, color)) {
                ++num_with_color;
            }
        }

        return num_with_color;
    }

    size_t numPredecessors(typename T::unitigmap_t const& unitig) override {
        size_t num_with_color = 0;
        for(auto& pred : unitig.getPredecessors()) {
            auto colorset = pred.getData()->getUnitigColors(pred);
            // Make mapping only include first or last k-mer (depending on whether we are on the reverse complement
            // or not)
            typename T::unitigmap_t last_kmer(pred);
            if(last_kmer.strand) {
                last_kmer.dist = last_kmer.len - 1;
            }

            last_kmer.len = 1;
            last_kmer.strand = true;

            if(colorset->contains(last_kmer, color)) {
                ++num_with_color;
            }
        }

        return num_with_color;
    }

private:
    size_t color;
};


/**
 * Special link annotator for reference genomes.
 *
 * This class *always* annotates the node in the graph containing the first k-mer of the sequence, such that it is
 * always possible to reconstruct the genome from it's starting node.
 *
 * @tparam T
 */
template<typename T>
class RefLinkAnnotator : public ColorAssociatedAnnotator<T> {
public:
    RefLinkAnnotator(T* _graph, LinkDB* _db, size_t _color) : ColorAssociatedAnnotator<T>(_graph, _db, _color) { }
    RefLinkAnnotator(RefLinkAnnotator const& o) = default;
    RefLinkAnnotator(RefLinkAnnotator&& o) noexcept = default;
    virtual ~RefLinkAnnotator() = default;

protected:
    bool nodeNeedsAnnotation(typename T::unitigmap_t const& unitig) override {
        if(!first_node_annotated) {
            first_node_annotated = true;
            return true;
        } else {
            return LinkAnnotator<T>::nodeNeedsAnnotation(unitig);
        }
    }

private:
    bool first_node_annotated = false;
};


/**
 * Read sequences from a FASTA file and add links for each sequence.
 */
template<typename T>
void addLinksFromFasta(LinkAnnotator<T>& annotator, vector<string> const& filepaths,
                       size_t batch_size=1e3) {
    atomic<bool> finished_reading(false);
    mutex queue_lock;
    queue<unique_ptr<vector<string>>> seq_blocks;
    condition_variable sequences_ready;

    thread link_creator_thread([&] () {
        while(true) {
            unique_lock<mutex> guard(queue_lock);
            sequences_ready.wait(guard, [&] () { return !seq_blocks.empty() || finished_reading; });

            if(seq_blocks.empty() && finished_reading) {
                cerr << "link_creator_thread :: finished.\n" << std::flush;
                return;
            }

            unique_ptr<vector<string>> sequences = std::move(seq_blocks.front());
            seq_blocks.pop();
            guard.unlock();

            for(string const& seq : *sequences) {
                annotator.addLinksFromSequence(seq);
            }

            cerr << "link_creator_thread :: processed block of " << sequences->size() << " sequences.\n" << flush;
        }
    });

    thread reader_thread([&] () {
        FileParser fp(filepaths);

        auto sequences = make_unique<vector<string>>();
        sequences->reserve(batch_size);
        string sequence;
        size_t file_ix = 0;
        size_t num_read = 0;
        while(fp.read(sequence, file_ix)) {
            ++num_read;
            sequences->emplace_back(sequence);

            if(num_read >= batch_size) {
                // When enough data read, push sequences on the queue. This way we don't have to lock the queue that
                // often.
                unique_lock<mutex> guard(queue_lock);
                cerr << "reader_thread :: read block of " << sequences->size() << " sequences.\n" << flush;
                seq_blocks.emplace(std::move(sequences));
                guard.unlock();

                sequences = make_unique<vector<string>>();
                sequences->reserve(batch_size);
                num_read = 0;
                sequences_ready.notify_one();
            }
        }

        // Push remaining sequences on the queue
        if(!sequences->empty()) {
            unique_lock<mutex> guard(queue_lock);
            cerr << "reader_thread :: read block of " << sequences->size() << " sequences.\n" << flush;
            seq_blocks.emplace(std::move(sequences));
            guard.unlock();

            sequences_ready.notify_one();
        }

        finished_reading = true;
        sequences_ready.notify_all();
    });

    link_creator_thread.join();
    reader_thread.join();
}

void define_MappingResult(py::module& m);

void define_LinkAnnotator(py::module& m);


}

#endif //PYFROST_LINKANNOTATOR_H
