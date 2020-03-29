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

    SECTION("This graph should have six unitigs") {
        int num_unitigs = 0;

        for (auto const& unitig : test_graph) {
            ++num_unitigs;
        }

        REQUIRE(num_unitigs == 6);
    }

    SECTION("Simplify shouldn't delete any short tips in this example") {
        test_graph.simplify();

        REQUIRE(test_graph.nbKmers() == 21);
    }

    SECTION("Check successors of an unitig") {
        auto start_unitig = test_graph.find(Kmer("ACTGA"), true).mappingToFullUnitig();

        auto unitig_str = start_unitig.referenceUnitigToString();
        if(!start_unitig.strand) {
            unitig_str = reverse_complement(unitig_str);
        }

        REQUIRE(unitig_str == "ACTGATTTCGA");

        size_t num_succ = 0;
        for(auto const& succ : start_unitig.getSuccessors()) {
            ++num_succ;

            auto succ_str = succ.referenceUnitigToString();
            if(!succ.strand) {
                succ_str = reverse_complement(succ_str);
            }

            REQUIRE(unitig_str.substr(unitig_str.size() - 4) == succ_str.substr(0, 4));
        }

        REQUIRE(num_succ == 2);

        // The reverse complement of the above unitig shouldn't have any successors
        auto start_unitig_rev = test_graph.find(Kmer("TCGAA"), true).mappingToFullUnitig();

        // This reverse complement unitig should be one of the successors of the original `start_unitig`.
        int matches = 0;
        for(auto const& succ : start_unitig.getSuccessors()) {
            if(succ == start_unitig_rev) {
                ++matches;
            }
        }

        REQUIRE(matches == 1);

        unitig_str = start_unitig_rev.referenceUnitigToString();
        if(!start_unitig_rev.strand) {
            unitig_str = reverse_complement(unitig_str);
        }

        REQUIRE(unitig_str == "TCGAAATCAGT");

        num_succ = 0;
        for(auto const& succ : start_unitig_rev.getSuccessors()) {
            ++num_succ;

            auto succ_str = succ.referenceUnitigToString();
            if(!succ.strand) {
                succ_str = reverse_complement(succ_str);
            }

            REQUIRE(unitig_str.substr(unitig_str.size() - 4) == succ_str.substr(0, 4));
        }

        REQUIRE(num_succ == 0);
    }
}

