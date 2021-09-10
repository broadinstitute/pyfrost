#ifndef PYFROST_LINKANNOTATOR_H
#define PYFROST_LINKANNOTATOR_H

#include "pyfrost.h"
#include "JunctionTree.h"
#include "LinkDB.h"

namespace pyfrost {

struct MappingResult {
    MappingResult() : mapping_start(0), mapping_end(0), num_junctions(0) { }
    MappingResult(MappingResult const& o) = default;
    MappingResult(MappingResult&& o) = default;

    /// Head k-mer of the first unitig this sequence maps to
    Kmer start_unitig;

    /// Head k-mer of the last unitig this sequence maps to
    Kmer end_unitig;

    /// First position in the sequence that has a k-mer that maps to the graph
    size_t mapping_start;

    /// Last position in the sequence that has a k-mer that maps to the graph
    size_t mapping_end;

    /// Number of junctions annotated thanks to this sequence
    size_t num_junctions;

    /// For each k-mer in the sequence, whether it was found in the graph or not
    std::vector<uint8_t> matches;

    /// For each unitig this sequence traversed, count how often we encountered it
    std::unordered_map<Kmer, size_t> unitig_visits;
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
    LinkAnnotator() = default;
    LinkAnnotator(LinkAnnotator const& o) = default;
    LinkAnnotator(LinkAnnotator&& o) noexcept = default;
    virtual ~LinkAnnotator() = default;

    /**
     * Thread the given sequence through the graph and annotate specific nodes with the junctions taken.
     *
     * @param graph The ccDBG
     * @param db Database to store link annotations
     * @param seq the sequence to thread through the graph
     */
    virtual MappingResult addLinksFromSequence(T& graph, LinkDB& db, std::string const& seq);

    /**
     * Thread the given sequence through the graph and annotate specific nodes with the junctions taken.
     *
     * If calling `addLinksFromSequence` multiple times, by default each sequence will be considered
     * independently, and previously annotated nodes will not be considered for the new sequence. In case of paired-end
     * reads, however, when processing the second read in a pair it might be worthwhile to continue adding links to
     * nodes annotated by the first read. In that case you can set `keep_nodes` to true.
     *
     * @param graph The ccDBG
     * @param db Database to store link annotations
     * @param seq the sequence to thread through the graph
     * @param keep_nodes Whether to keep the nodes to annotate from a previous call
     */
    virtual MappingResult addLinksFromSequence(T& graph, LinkDB& db, std::string const& seq, bool keep_nodes);


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

    std::vector<std::shared_ptr<JunctionTreeNode>> nodes_to_annotate;
};

template<typename T>
bool LinkAnnotator<T>::nodeNeedsAnnotation(typename T::unitigmap_t const& unitig) {
    for(auto& succ : unitig.getSuccessors()) {
        if(succ.getPredecessors().cardinality() > 1) {
            return true;
        }
    }

    for(auto& pred : unitig.getPredecessors()) {
        if(pred.getSuccessors().cardinality() > 1) {
            return true;
        }
    }

    return false;
}

template<typename T>
MappingResult LinkAnnotator<T>::addLinksFromSequence(T& graph, LinkDB& db, string const& seq) {
    return LinkAnnotator<T>::addLinksFromSequence(graph, db, seq, false);
}

template<typename T>
MappingResult LinkAnnotator<T>::addLinksFromSequence(T& graph, LinkDB& db, string const& seq, bool keep_nodes) {
    if(!keep_nodes) {
        resetNodesToAnnotate();
    }

    bool first_unitig_found = false;
    MappingResult mapping;
    mapping.start_unitig.set_empty();
    mapping.end_unitig.set_empty();
    mapping.matches.resize(seq.length() - Kmer::k + 1, 0);

    KmerIterator kmer_iter(seq.c_str()), kmer_end;
    for(; kmer_iter != kmer_end; ++kmer_iter) {
        Kmer kmer;
        size_t pos;
        tie(kmer, pos) = *kmer_iter;

        //cerr << "\nKmer " << kmer.toString() << ", pos: " << pos << endl;

        // Check if we're at the end of an unitig
        auto umap = graph.find(kmer);
        if(umap.isEmpty) {
            // Kmer doesn't exist in the graph, move to next
            continue;
        }

        auto unitig = umap.mappingToFullUnitig();
        auto unitig_kmer = unitig.getMappedHead();
        //cerr << "- on unitig: " << unitig_kmer.toString() << endl;

        if(!first_unitig_found) {
            // First k-mer that is present in the graph
            mapping.start_unitig = unitig_kmer;
            mapping.mapping_start = pos;
            first_unitig_found = true;
        }

        mapping.matches[pos] = true;
        mapping.mapping_end = pos;
        mapping.end_unitig = unitig_kmer;

        auto it = mapping.unitig_visits.find(unitig_kmer);
        if(it != mapping.unitig_visits.end()) {
            ++(it->second);
        } else {
            mapping.unitig_visits[unitig_kmer] = 1;
        }

        if(nodeNeedsAnnotation(unitig)) {
            nodes_to_annotate.push_back(db.createOrGetTree(unitig.getMappedTail()));
        }

        // Move to the end of the unitig, and by definition we will not encounter any branch points.
        size_t unitig_len = unitig.size - graph.getK() + 1; // unitig length in num k-mers
        size_t oriented_pos = kmerPosOriented(umap);
        size_t diff_to_unitig_end = unitig_len - oriented_pos - 1;

//        cerr << "- len: " << unitig_len << ", pos: " << umap.dist << ", oriented pos: " << oriented_pos
//             << ", diff: " << diff_to_unitig_end << endl;

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
//            cerr << "End of sequence." << endl;
            break;
        }

        if(kmer != unitig.getMappedTail()) {
            // The k-mer doesn't match the unitig tail anymore. The given sequence (a read) likely contained an error
            // that was removed from the graph, and now there's no unambiguous path anymore, so we quit.
//            cerr << "Kmer doesn't match unitig tail!" << kmer.toString() << " " << unitig.getMappedTail().toString()
//                 << endl;
            break;
        }

        size_t edge_pos = pos + Kmer::k;
        if(edge_pos >= seq.length()) {
            // At the end of the sequence, so no edge to inspect
//            cerr << "No sequence left anymore" << endl;
            break;
        }

//        cerr << "- at unitig end: " << unitig.getMappedTail().toString() << endl;

        if(unitig.getSuccessors().cardinality() > 1) {
            char edge = toupper(seq[edge_pos]);
            if(!(edge == 'A' || edge == 'C' || edge == 'G' || edge == 'T')) {
                // Invalid sequence, quit
                break;
            }

            Kmer succ_kmer = kmer.forwardBase(edge);
            bool found_succ = false;

            // Let's see if there's a successor that matches with the sequence
            for(auto& succ : umap.getSuccessors()) {
                if(succ.getMappedHead() == succ_kmer) {
                    // Add edge choice to each tree
                    std::transform(
                        nodes_to_annotate.begin(), nodes_to_annotate.end(), nodes_to_annotate.begin(),
                        [edge] (std::shared_ptr<JunctionTreeNode> node) {
                            return node->addEdge(edge);
                        });
                    ++mapping.num_junctions;

                    found_succ = true;
//                    cerr << "- " << edge << ", successor: " << succ_kmer.toString() << ", nodes annotated: "
//                         << nodes_to_annotate.size() << endl;
                    break;
                }
            }

            // Couldn't find a valid successor in the graph, abort
            if(!found_succ) {
                break;
            }
        }
    }

    return mapping;
}

/**
 * Special link annotator for reference genomes.
 *
 * This class *always* annotates the node in the graph containing the first k-mer of the sequence, such that it is
 * always possible to reconstruct the genome from it's starting node.
 *
 * @tparam T
 */
template<typename T>
class RefLinkAnnotator : public LinkAnnotator<T> {
public:
    RefLinkAnnotator() = default;
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
void addLinksFromFasta(T& graph, LinkDB& db, LinkAnnotator<T>& annotator, vector<string> const& filepaths,
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
                annotator.addLinksFromSequence(graph, db, seq);
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
