/**
 * A pretty simple k-mer counter.
 *
 * Based on ideas of Heng Li:
 * https://github.com/lh3/kmer-cnt
 *
 * This k-mer counter uses an ensemble of hash tables, and the k-mer minimizer hash
 * is used as index which hash table to use. This is especially useful when multithreading
 * the counting of k-mers.
 *
 * We created this k-mer counter because Bifrost doesn't provide k-mer counts on its own
 * and well engineered existing counters were all licensed under GPL-3.
 */

#ifndef PYFROST_KMERCOUNTER_H
#define PYFROST_KMERCOUNTER_H

#include "Pyfrost.h"

#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <robin_hood.h>
#include <File_Parser.hpp>

#include "Kmer.h"

namespace pyfrost {

class KmerCounter {
public:
    KmerCounter(size_t k, size_t g=0, bool canonical=true, size_t num_threads=2, size_t table_bits=10,
        size_t read_block_size=(1 << 20));

    KmerCounter& countKmersFiles(std::vector<std::string> const& files);
    KmerCounter& countKmers(std::string const& str);

    uint16_t query(Kmer const& kmer);

    uint64_t getNumKmers() const {
        return num_kmers.load();
    }

    uint64_t getUniqueKmers() const {
        return num_unique.load();
    }

private:
    void counterThread();
    void setKmerGmer();

    size_t k;
    size_t g;
    bool canonical;
    size_t num_threads;
    size_t read_block_size;

    std::vector<robin_hood::unordered_map<Kmer, uint16_t>> tables;
    std::vector<std::mutex> table_locks;
    std::atomic<uint64_t> num_kmers;
    std::atomic<uint64_t> num_unique;

    std::atomic<bool> finished_reading;
    std::mutex queue_lock;
    std::condition_variable sequence_ready;
    std::queue<std::string> seq_queue;

};

void define_KmerCounter(py::module& m);

}





#endif //PYFROST_KMERCOUNTER_H
