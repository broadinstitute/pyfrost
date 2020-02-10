#include <iostream>
#include <catch2/catch.hpp>
#include <CompactedDBG.hpp>

TEST_CASE( "Test simple bubble detection", "[bubbles]") {
    CompactedDBG<> test_graph = CompactedDBG<>();

    std::string seq1 = "GTCTAGCATGCATCGATGCTAGCTAGCTAGCTA";
    std::string seq2 =

    bool result = test_graph.add(seq1);
    REQUIRE(result);

    bool result2 = test_graph.add("GTCTAGCATGCATCGATGCTAGCTAGCTAGCTA");
    REQUIRE(result2);
}
