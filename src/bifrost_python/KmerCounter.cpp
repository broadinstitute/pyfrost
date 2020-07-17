#include <thread>
#include <cereal/cereal.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

#include "KmerCounter.h"

using std::vector;
using std::string;
using std::thread;
using std::lock_guard;
using std::unique_lock;
using std::mutex;


namespace pyfrost {

KmerCounter::KmerCounter(size_t _k, size_t _g, bool _canonical, size_t _num_threads, size_t _table_bits,
    size_t _read_block_size) :
    k(_k), g(_g), canonical(_canonical), num_threads(_num_threads), read_block_size(_read_block_size),
    tables(1 << _table_bits), table_locks(1 << _table_bits),
    num_kmers(0), num_unique(0), finished_reading(false)
{
    setKmerGmer();
}

KmerCounter::KmerCounter(KmerCounter const& o) :
    k(o.k), g(o.g), canonical(o.canonical), num_threads(o.num_threads), read_block_size(o.read_block_size),
    tables(o.tables), table_locks(o.tables.size()), num_kmers(o.num_kmers.load()), num_unique(o.num_unique.load()),
    finished_reading(o.finished_reading.load())
{
    setKmerGmer();
}

KmerCounter::KmerCounter(KmerCounter&& o) :
    k(o.k), g(o.g), canonical(o.canonical), num_threads(o.num_threads), read_block_size(o.read_block_size),
    tables(std::move(o.tables)), table_locks(o.tables.size()), num_kmers(o.num_kmers.load()),
    num_unique(o.num_unique.load()), finished_reading(o.finished_reading.load())
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
    unique_lock<mutex> guard(queue_lock);
    seq_queue.emplace(str);
    guard.unlock();
    finished_reading = true;

    thread counter_thread(&KmerCounter::counterThread, this);
    sequence_ready.notify_one();

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

                unique_lock<mutex> table_guard(table_locks[table_ix]);
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
                table_guard.unlock();
            }
        }
    }
}

uint16_t KmerCounter::query(const Kmer &qry) const
{
    Kmer kmer = canonical ? qry.rep() : qry;
    string kmer_str = kmer.toString();

    minHashKmer<RepHash> it_min(kmer_str.c_str(), Kmer::k, Minimizer::g, RepHash(), true), it_min_end;
    for(; it_min != it_min_end; ++it_min) {
        uint64_t minimizer_hash = it_min.getHash();
        size_t table_ix = minimizer_hash % tables.size();

        auto it = tables[table_ix].find(kmer);
        if(it != tables[table_ix].end()) {
            return it->second;
        }
    }

    return 0;
}

uint16_t KmerCounter::query(char const* qry) const {
    return query(Kmer(qry));
}

template<typename Archive>
void KmerCounter::save(Archive& ar) const {
    ar(k, g, canonical, tables, num_kmers.load(), num_unique.load());
}

template<typename Archive>
void KmerCounter::load(Archive& ar) {
    ar(k, g, canonical);
    setKmerGmer();

    size_t _num_kmers = 0;
    size_t _num_unique = 0;
    ar(tables, _num_kmers, _num_unique);

    num_kmers = _num_kmers;
    num_unique = _num_unique;

    table_locks = vector<mutex>(tables.size());
}


void define_KmerCounter(py::module& m) {
    auto py_KmerCounter = py::class_<KmerCounter>(m, "KmerCounter")
        .def(py::init<size_t, size_t, bool, size_t, size_t, size_t>(),
            py::arg("k"), py::arg("g") = 0, py::arg("canonical") = true,
            py::arg("num_threads") = 2, py::arg("table_bits") = 10, py::arg("read_block_size") = (1 << 20))
        .def("count_kmers", &KmerCounter::countKmers)
        .def("count_kmers_files", &KmerCounter::countKmersFiles)
        .def("query", py::overload_cast<Kmer const&>(&KmerCounter::query, py::const_))
        .def("query", py::overload_cast<char const*>(&KmerCounter::query, py::const_))

        .def("__getitem__", py::overload_cast<Kmer const&>(&KmerCounter::query, py::const_), py::is_operator())
        .def("__getitem__", py::overload_cast<char const*>(&KmerCounter::query, py::const_), py::is_operator())

        .def("__iter__", [] (KmerCounter& self) {
            return py::make_key_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        .def("items", [] (KmerCounter& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        .def("__len__", [] (KmerCounter const& self) {
            return self.getUniqueKmers();
        })
        .def_property_readonly("num_kmers", &KmerCounter::getNumKmers)
        .def_property_readonly("num_unique", &KmerCounter::getUniqueKmers)

        .def("save", [] (KmerCounter& self, string const& filepath) {
            std::ofstream ofile(filepath);
            cereal::BinaryOutputArchive archive(ofile);

            archive(cereal::make_nvp("kmercounter", self));
        })

        .def_static("from_file", [] (string const& filepath) {
            KmerCounter kmer_counter;

            std::ifstream ifile(filepath);
            cereal::BinaryInputArchive archive(ifile);

            archive(kmer_counter);

            return kmer_counter;
        });

    auto Mapping = py::module::import("collections.abc").attr("Mapping");
    auto Set = py::module::import("collections.abc").attr("Set");
    py_KmerCounter.attr("__bases__") = py::make_tuple(Mapping, Set).attr("__add__")(py_KmerCounter.attr("__bases__"));
}

}


//
// Cereal support
// --------------
// Make sure we can serialize Kmer objects and the robin_hood unordered_map
//

template<typename Archive>
std::string save_minimal(Archive& ar, Kmer const& kmer)
{
    stringstream bytes;
    kmer.write(bytes);

    return bytes.str();
}

template<typename Archive>
void load_minimal(Archive& ar, Kmer& kmer, std::string const& value)
{
    istringstream stream(value);
    kmer.read(stream);
}

template<typename Archive, typename K, typename V, typename C, typename A>
void save(Archive& ar, robin_hood::unordered_map<K, V, C, A> const& map)
{
    ar(map.size());
    for(auto const& e : map) {
        ar(e.first, e.second);
    }
}

template<typename Archive, typename K, typename V, typename C, typename A>
void load(Archive& ar, robin_hood::unordered_map<K, V, C, A>& map)
{
    map.clear();

    size_t num_entries;
    ar(num_entries);
    map.reserve(num_entries);

    for(int i = 0; i < num_entries; ++i) {
        K key;
        V value;
        ar(key);
        ar(value);

        map.emplace(std::move(key), std::move(value));
    }
}
