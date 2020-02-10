#include <CompactedDBG.hpp>

#include <iostream>

int main(int argc, char** argv) {
    CompactedDBG<> test = CompactedDBG<>();

    test.add("AAAAAAAAAAAAAAAAAGTCTAGCATGCATCGATGCTAGCTAGCTAGCTA");
    test.add("AAAAAAAAAAAAAAAAAGTCTAGCATGCATCGAGGCTAGCTAGCTAGCTA");

    for(auto& unitig : test) {
        std::cout << "Unitig: " << unitig.referenceUnitigToString() << std::endl;
    }

    return 0;
}