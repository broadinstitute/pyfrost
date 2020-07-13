#include <thread>
#include "KmerCounter.h"

using std::vector;
using std::string;
using std::thread;
using std::lock_guard;
using std::unique_lock;
using std::mutex;

namespace std {

template<>
class hash<Kmer> {
public:
    size_t operator()(Kmer const& kmer) const {
        return kmer.hash();
    }
};

}

namespace pyfrost {

KmerCounter::KmerCounter(size_t _k, size_t _g, bool _canonical, size_t _num_threads, size_t _table_bits,
    size_t _read_block_size) :
    k(_k), g(_g), canonical(_canonical), num_threads(_num_threads), read_block_size(_read_block_size),
    tables(1 << _table_bits), table_locks(1 << _table_bits),
    num_kmers(0), num_unique(0), finished_reading(false)
{
    setKmerGmer();
}

void KmerCounter::setKmerGmer() {
    Kmer::set_k(k);

    if (g == 0) {
        if(k >= 15) {
            g = k - DEFAULT_G_DEC1;
        } else if(k >= 7) {
            g = k - DEFAULT_G_DEC2;
        } else {
            g = k - 2;
        }
    }

    Minimizer::set_g(g);
}

KmerCounter& KmerCounter::countKmers(std::string const& str)
{
    finished_reading = false;
    thread counter_thread(&KmerCounter::counterThread, this);

    unique_lock<mutex> guard(queue_lock);
    seq_queue.emplace(str);
    guard.unlock();

    sequence_ready.notify_one();

    finished_reading = true;
    sequence_ready.notify_all();

    counter_thread.join();

    return *this;
}

KmerCounter& KmerCounter::countKmersFiles(std::vector<std::string> const& files)
{
    vector<thread> counter_threads;
    finished_reading = false;

    for(int i = 0; i < num_threads; ++i) {
        counter_threads.emplace_back(&KmerCounter::counterThread, this);
    }

    std::thread reader_thread([&] () {
        FileParser fp(files);

        std::vector<string> sequences;
        string sequence;
        size_t file_ix = 0;
        size_t nucleotides_read = 0;
        while(fp.read(sequence, file_ix)) {
            nucleotides_read += sequence.size();
            sequences.emplace_back(sequence);

            if(nucleotides_read >= read_block_size) {
                // When enough data read, push sequences on the queue. This way we don't have to lock the queue that
                // often.
                unique_lock<mutex> guard(queue_lock);
                for(auto& seq : sequences) {
                    seq_queue.emplace(std::move(seq));
                }
                guard.unlock();

                sequences.clear();
                nucleotides_read = 0;
                sequence_ready.notify_one();
            }
        }

        // Push remaining sequences on the queue
        if(!sequences.empty()) {
            unique_lock<mutex> guard(queue_lock);
            for(auto& seq : sequences) {
                seq_queue.emplace(std::move(seq));
            }
            guard.unlock();

            sequences.clear();
            nucleotides_read = 0;
            sequence_ready.notify_one();

        }

        finished_reading = true;
        sequence_ready.notify_all();
    });

    for(auto& t : counter_threads) {
        t.join();
    }

    reader_thread.join();

    return *this;
}


void KmerCounter::counterThread()
{
    while(true) {
        // Wait until sequences are available in the queue
        unique_lock<mutex> guard(queue_lock);
        sequence_ready.wait(guard, [&] () { return !seq_queue.empty() || finished_reading; });

        if(seq_queue.empty() && finished_reading) {
            return;
        }

        // Obtain block of sequences to count
        vector<string> sequences;
        size_t nucleotides_read = 0;
        while(nucleotides_read < read_block_size && !seq_queue.empty()) {
            sequences.emplace_back(std::move(seq_queue.front()));
            seq_queue.pop();
        }
        guard.unlock();

        for(auto& sequence : sequences) {
            KmerIterator kmer_iter(sequence.c_str()), kmer_end;
            minHashIterator<RepHash> it_min(sequence.c_str(), sequence.size(), Kmer::k, Minimizer::g, RepHash(), true);

            for(; kmer_iter != kmer_end; ++kmer_iter) {
                std::pair<Kmer, int> const p = *kmer_iter; // K-mer hash and position in sequence
                Kmer kmer = canonical ? p.first.rep() : p.first;
                ++num_kmers;

                // Move minimizer position, also takes into account if one or more k-mer were jumped because contained
                // non-ACGT char.
                // Minimizer hash serves as hash table index
                it_min += (p.second - it_min.getKmerPosition());
                uint64_t minimizer_hash = it_min.getHash();
                size_t table_ix = minimizer_hash % tables.size();

                lock_guard<mutex> table_guard(table_locks[table_ix]);
                uint16_t curr_count = 0;
                if(tables[table_ix].contains(kmer)) {
                    curr_count = tables[table_ix][kmer];
                } else {
                    ++num_unique;
                }

                // Check if counter is not saturated
                if(curr_count < 0xFFFF){
                    ++curr_count;
                }

                tables[table_ix].insert_or_assign(kmer, curr_count);
            }
        }
    }
}

uint16_t KmerCounter::query(const Kmer &qry)
{
    Kmer kmer = canonical ? qry.rep() : qry;
    string kmer_str = kmer.toString();

    minHashKmer<RepHash> it_min(kmer_str.c_str(), Kmer::k, Minimizer::g, RepHash(), true), it_min_end;
    for(; it_min != it_min_end; ++it_min) {
        uint64_t minimizer_hash = it_min.getHash();
        size_t table_ix = minimizer_hash % tables.size();

        if(tables[table_ix].contains(kmer)) {
            return tables[table_ix][kmer];
        }
    }

    return 0;
}





void define_KmerCounter(py::module& m) {
    py::class_<KmerCounter>(m, "KmerCounter")
        .def(py::init<size_t, size_t, bool, size_t, size_t, size_t>(),
            py::arg("k"), py::arg("g") = 0, py::arg("canonical") = true,
            py::arg("num_threads") = 2, py::arg("table_bits") = 10, py::arg("read_block_size") = (1 << 20))
        .def("count_kmers", &KmerCounter::countKmers)
        .def("count_kmers_files", &KmerCounter::countKmersFiles)
        .def("query", &KmerCounter::query)

        .def_property_readonly("num_kmers", &KmerCounter::getNumKmers)
        .def_property_readonly("num_unique", &KmerCounter::getUniqueKmers);
}

}


