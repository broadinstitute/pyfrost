#include <iostream>
#include <unordered_set>

#include <catch2/catch.hpp>
#include <CompactedDBG.hpp>
#include <Kmer.hpp>

TEST_CASE("Test McCortex example", "[dbg_construction]") {
    CDBG_Build_opt opt{};

    opt.k = 5;
    opt.g = 3;
    opt.filename_ref_in.emplace_back("data/mccortex.fasta");

    CompactedDBG<> test_graph = CompactedDBG<>(opt.k, opt.g);
    test_graph.build(opt);

    REQUIRE(test_graph.nbKmers() == 21);

    SECTION("This graph should have four unitigs") {
        int num_unitigs = 0;

        for (auto const& unitig : test_graph) {
            ++num_unitigs;

            std::cerr << "Unitig: " << unitig.referenceUnitigToString() << std::endl;
            std::cerr << "Rev compl: " << reverse_complement(unitig.referenceUnitigToString()) << std::endl;

            for (auto &p : unitig.getPredecessors()) {
                std::cerr << "- Predecessor: " << p.referenceUnitigToString() << std::endl;
            }

            for (auto &s : unitig.getSuccessors()) {
                std::cerr << "- Successor: " << s.referenceUnitigToString() << std::endl;
            }
        }

        REQUIRE(num_unitigs == 4);
    }

    SECTION("Simplify shouldn't delete any short tips in this example") {
        test_graph.simplify();

        REQUIRE(test_graph.nbKmers() == 21);
    }
}

TEST_CASE( "Test simple bubble detection", "[bubbles]") {
    CDBG_Build_opt opt{};

    opt.k = 5;
    opt.g = 3;
    opt.filename_ref_in.emplace_back("data/bubble_test.fasta");

    CompactedDBG<> test_graph = CompactedDBG<>(opt.k, opt.g);
    test_graph.build(opt);

    for(auto const& unitig : test_graph) {
        std::cerr << "Unitig: " << unitig.referenceUnitigToString() << std::endl;

        for(auto& p : unitig.getPredecessors()) {
            std::cerr << "Predecessor: " << p.referenceUnitigToString() << std::endl;
        }

        for(auto& s : unitig.getSuccessors()) {
            std::cerr << "Successor: " << s.referenceUnitigToString() << std::endl;
        }
    }

    Kmer start = Kmer("TAATG");
    auto start_unitig = test_graph.find(start, true);

    std::cerr << "Start unitig: " << start_unitig.referenceUnitigToString() << std::endl;

    test_graph.write("test_bubble");
}
