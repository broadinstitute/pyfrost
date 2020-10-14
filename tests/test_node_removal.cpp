#include <iostream>
#include <unordered_set>

#include <catch2/catch.hpp>
#include <ColoredCDBG.hpp>
#include <Kmer.hpp>

TEST_CASE("Test node/color removal", "[node_color_removal]") {
    CCDBG_Build_opt opt;
    opt.filename_graph_in = "data/F11-frags.gfa";
    opt.filename_colors_in = "data/F11-frags.bfg_colors";

    ColoredCDBG<> ccdbg(opt.k, opt.g);
    ccdbg.read("data/F11-frags.gfa", "data/F11-frags.bfg_colors", 2);

    for(auto const& um : ccdbg) {
        auto colorset = um.getData()->getUnitigColors(um);
        REQUIRE(colorset != nullptr);

        int num_colors = 0;
        for(auto it = colorset->begin(um); it != colorset->end(); ++it) {
            ++num_colors;
        }

        REQUIRE(num_colors > 0);
    }

    ifstream to_remove;
    to_remove.open("data/to_remove.txt");

    const int SAMPLE_COLOR_ID = 0;

    std::vector<Kmer> unitigs_to_remove;
    while(to_remove.good()) {
        std::string node;
        to_remove >> node;

        Kmer kmer = Kmer(node.c_str());
        auto um = ccdbg.find(kmer, true).mappingToFullUnitig();
        REQUIRE(!um.isEmpty);

        um.strand = true;
        um.dist = 0;

        auto colorset = um.getData()->getUnitigColors(um);
        colorset->remove(um, SAMPLE_COLOR_ID);

        if(colorset->size(um) == 0) {
            unitigs_to_remove.emplace_back(um.getUnitigHead());
        }
    }

    for(auto const& um : ccdbg) {
        auto colorset = um.getData()->getUnitigColors(um);
        REQUIRE(colorset != nullptr);
    }

    for(auto const& kmer : unitigs_to_remove) {
        auto um = ccdbg.find(kmer, true).mappingToFullUnitig();
        if(um.isEmpty) {
            // Might be already removed because of rev compl k-mer
            continue;
        }

        um.strand = true;
        ccdbg.remove(um);
    }

    int num_unitigs_without_colors = 0;
    for(auto const& um : ccdbg) {
        auto colorset = um.getData()->getUnitigColors(um);
        REQUIRE(colorset != nullptr);

        int num_colors = 0;
        for(auto it = colorset->begin(um); it != colorset->end(); ++it) {
            ++num_colors;
        }

        if(num_colors == 0) {
            ++num_unitigs_without_colors;
        }
    }

    REQUIRE(num_unitigs_without_colors == 0);

    ccdbg.write("cleaned");

    ColoredCDBG<> ccdbg2(opt.k, opt.g);
    ccdbg2.read("cleaned.gfa", "cleaned.bfg_colors", 2);

    for(auto const& um : ccdbg) {
        auto colorset = um.getData()->getUnitigColors(um);
        REQUIRE(colorset != nullptr);

        int num_colors = 0;
        for(auto it = colorset->begin(um); it != colorset->end(); ++it) {
            ++num_colors;
        }

        REQUIRE(num_colors > 0);
    }
}

