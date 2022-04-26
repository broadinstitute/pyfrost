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

#include "pyfrost.h"
#include "Kmer.h"
#include "Serialize.h"

#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <robin_hood.h>
#include <File_Parser.hpp>

namespace pyfrost {

using std::pair;
using std::vector;

class KmerCounter;

typedef uint16_t kmercount_t;
using KmerCountMap = robin_hood::unordered_map<Kmer, kmercount_t>;

class KmerCounterIterator {
public:
    using iteratory_category = std::input_iterator_tag;
    using value_type = pair<Kmer, kmercount_t>;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using pointer = value_type*;

    KmerCounterIterator(vector<KmerCountMap>::iterator const& _table_iter,
        vector<KmerCountMap>::iterator const& _table_end) :
        curr_table(_table_iter), table_end(_table_end)
    {
        if(curr_table != table_end) {
            do {
                curr_kmer = curr_table->begin();
                kmer_end = curr_table->end();

                if(curr_kmer == kmer_end) {
                    ++curr_table;
                }
            } while(curr_kmer == kmer_end && curr_table != table_end);
        }
    }

    KmerCounterIterator(KmerCounterIterator const& o) = default;
    KmerCounterIterator(KmerCounterIterator&& o) = default;

    KmerCounterIterator& operator++() {
        if(curr_table == table_end) {
            return *this;
        }

        ++curr_kmer;
        if(curr_kmer == kmer_end) {
            ++curr_table;
            if(curr_table != table_end) {
                do {
                    curr_kmer = curr_table->begin();
                    kmer_end = curr_table->end();

                    if(curr_kmer == kmer_end) {
                        ++curr_table;
                    }
                } while(curr_kmer == kmer_end && curr_table != table_end);
            }
        }

        return *this;
    }

    KmerCounterIterator operator++(int) {
        KmerCounterIterator tmp(*this);

        operator++();
        return tmp;
    }

    value_type operator*() {
        if(curr_table == table_end) {
            Kmer tmp;
            tmp.set_empty();

            return {tmp, 0u};
        }

        // Make std::pair to ensure automatic conversion to a py::tuple
        return make_pair(curr_kmer->first, curr_kmer->second);
    }

    robin_hood::pair<const Kmer, kmercount_t>* operator->() {
        if(curr_table == table_end) {
            return nullptr;
        }

        return curr_kmer.operator->();
    }

    bool operator==(KmerCounterIterator const& o) {
        if(curr_table == table_end && o.curr_table == o.table_end) {
            return true;
        } else {
            return curr_table == o.curr_table && curr_kmer == o.curr_kmer;
        }
    }

    bool operator!=(KmerCounterIterator const& o) {
        return !operator==(o);
    }

private:
    vector<KmerCountMap>::iterator curr_table;
    vector<KmerCountMap>::iterator table_end;
    KmerCountMap::iterator curr_kmer;
    KmerCountMap::iterator kmer_end;

};

class KmerCounter {
public:
    KmerCounter(size_t k=DEFAULT_K, size_t g=0, bool canonical=true, size_t num_threads=2, size_t table_bits=10,
        size_t batch_size=100000);

    KmerCounter(KmerCounter const& o);
    KmerCounter(KmerCounter&& o);

    KmerCounter& countKmersFiles(std::vector<std::string> const& files);
    KmerCounter& countKmers(std::string const& str);

    kmercount_t query(char const* qry) const;
    kmercount_t query(Kmer const& qry) const;

    uint64_t getNumKmers() const {
        return num_kmers.load();
    }

    uint64_t getUniqueKmers() const {
        return num_unique.load();
    }

    kmercount_t getMaxCount() const {
        return max_count.load();
    }

    std::vector<uint64_t> getFrequencySpectrum();

    template<typename Archive>
    void save(Archive& ar) const {
        ar(k, g, canonical, tables, num_kmers.load(), num_unique.load(), max_count.load());
    }

    template<typename Archive>
    void load(Archive& ar) {
        ar(k, g, canonical);
        setKmerGmer();

        size_t _num_kmers = 0;
        size_t _num_unique = 0;
        kmercount_t _max_count = 0;
        ar(tables, _num_kmers, _num_unique, _max_count);

        num_kmers = _num_kmers;
        num_unique = _num_unique;
        max_count = _max_count;

        table_locks = vector<mutex>(tables.size());
    }

    KmerCounterIterator begin() {
        return {tables.begin(), tables.end()};
    }

    KmerCounterIterator end() {
        return {tables.end(), tables.end()};
    }

private:
    void counterThread();
    void setKmerGmer();

    size_t k;
    size_t g;
    bool canonical;
    size_t num_threads;
    size_t batch_size;

    std::vector<KmerCountMap> tables;
    std::vector<std::mutex> table_locks;
    std::atomic<uint64_t> num_kmers;
    std::atomic<uint64_t> num_unique;
    std::atomic<kmercount_t> max_count;

    std::atomic<bool> finished_reading;
    std::mutex queue_lock;
    std::condition_variable sequence_ready;
    std::queue<std::unique_ptr<std::vector<std::string>>> seq_queue;

};

void define_KmerCounter(py::module& m);

}

#endif //PYFROST_KMERCOUNTER_H
