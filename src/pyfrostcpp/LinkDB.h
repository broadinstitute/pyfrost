#ifndef PYFROST_LINKDB_H
#define PYFROST_LINKDB_H

#include <string>
#include <thread>
#include <queue>
#include <robin_hood.h>

#include "pyfrost.h"
#include "JunctionTree.h"

namespace pyfrost {

using junction_tree_map = robin_hood::unordered_map<Kmer, std::shared_ptr<JunctionTreeNode>>;

class LinkDB {
public:
    LinkDB() {}

    LinkDB(LinkDB const &o) = default;
    LinkDB(LinkDB &&o) noexcept = default;

    virtual ~LinkDB() = default;

    LinkDB& operator=(LinkDB const& o) = default;

    /**
     * Check if the database has links for a given kmer
     */
    virtual bool hasLinks(Kmer const& kmer) const = 0;

    /**
     * Get a JunctionTreeNode representing all links for a given k-mer
     */
    virtual std::shared_ptr<JunctionTreeNode> getLinks(Kmer const& kmer) = 0;

    /**
     * Total number of JunctionTrees in this database
     */
    virtual size_t numTrees() const = 0;

    /**
     * Return all <kmer, JunctionTreeNode> instances in this database.
     */
    virtual junction_tree_map &getJunctionTrees() = 0;

    /**
     * Create a new junction tree for the given k-mer in the database, or return the existing one if available.
     */
    virtual std::shared_ptr<JunctionTreeNode> createOrGetTree(Kmer const& kmer) = 0;

};

void define_LinkDB(py::module& m);

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
           : unitig.dist > 0 ? (unitig.size - unitig.getGraph()->getK() - unitig.dist) : 0;
}

/**
 * Thread a given `sequence` through the graph and annotate visited k-mers before junctions with a junction tree, to be
 * stored in the given `LinkDB` object.
 *
 * @tparam T Type of graph (CompactedDBG, ColoredCDBG, ...)
 * @param graph The graph object
 * @param db LinkDB object used to store junction trees.
 * @param sequence DNA sequence to thread through the graph
 */
template<typename T>
void addLinksFromSequence(T& graph, LinkDB& db, const std::string& sequence)
{
    KmerIterator kmer_iter(sequence.c_str()), kmer_end;
    // Keep track through which branch points this sequence is threaded through, to easily update previous junction
    // trees as the threading process progresses.
    auto nodes_to_annotate = std::vector<std::shared_ptr<JunctionTreeNode>>();

    typename T::unitigmap_t prev_unitig;
    for(; kmer_iter != kmer_end; ++kmer_iter) {
        std::pair<Kmer, int> p = *kmer_iter;
        auto umap = graph.find(p.first);

        // Check if k-mer exists in graph
        if(umap.isEmpty) {
            continue;
        }

        auto unitig = umap.mappingToFullUnitig();
        // If any of its successors have an in-degree > 1, then this node should be annotated with a junction tree.
        for(auto& succ : unitig.getSuccessors()) {
            if(succ.getPredecessors().cardinality() > 1) {
                nodes_to_annotate.push_back(db.createOrGetTree(unitig.getMappedTail()));
                break;
            }
        }

        prev_unitig = unitig;
        Kmer unitig_head = unitig.getMappedHead();

        size_t unitig_len = umap.size - graph.getK() + 1;
        size_t oriented_pos = kmerPosOriented(umap);
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
            break;
        }

        // We should be at the end of a unitig, so either we have a branch point, or our successor has a in-degree > 1
        // Make sure we have enough sequence left, and see what path the sequence takes through the graph
        size_t next_pos = p.second + Kmer::k;
        if(next_pos >= sequence.length()) {
            // We're already at the end of the link sequence, quit
            break;
        }

        if(unitig.getSuccessors().cardinality() > 1) {
            // We have found an unitig with two or more outgoing edges, let's see which branch
            // the sequence takes
            char edge = std::toupper(sequence[next_pos]);
            if(!(edge == 'A' || edge == 'C' || edge == 'G' || edge == 'T')) {
                // Invalid character in sequence, continue, KmerIterator will take care of skipping ahead
                continue;
            }

            Kmer next_kmer = p.first.forwardBase(edge);
            bool found_succ = false;

            // Let's see if there's an successor that matches with the link sequence
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

/**
 * Read sequences from a file and add them as links
 */
template<typename T>
void addLinksFromFile(T& graph, LinkDB& db, vector<string> const& filepaths, size_t batch_size=1e5) {
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
                addLinksFromSequence(graph, db, seq);
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

}

#endif //PYFROST_LINKDB_H
