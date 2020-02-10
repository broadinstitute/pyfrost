#include <catch2/catch.hpp>
#include <CompactedDBG.hpp>

TEST_CASE( "Test simple bubble detection", "[bubbles]") {
    CompactedDBG<> test_graph = CompactedDBG<>();

    bool result = test_graph.add("GTCTAGCATGCATCGATGCTAGCTAGCTAGCTA");
    REQUIRE(result);

    bool result2 = test_graph.add("GTCTAGCATGCATCGATGCTAGCTAGCTAGCTA");
    REQUIRE(result2);
}
