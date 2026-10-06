// dinara-align from C++: the edit distance and an optimal alignment's CIGAR for two sequences.
//
//     pixi run build-c
//     c++ -std=c++17 -I build/c c/example.cpp -L build/c -ldinara -Wl,-rpath,build/c -o example
//     ./example ACGTACGTTTGCA ACGTCGTTTTGCA
#include <iostream>

#include "dinara.h"

int main(int argc, char **argv) {
    std::string_view first = argc > 2 ? argv[1] : "ACGTACGTTTGCA";
    std::string_view second = argc > 2 ? argv[2] : "ACGTCGTTTTGCA";
    try {
        dinara::Alignment aligned = dinara::edit_cigar(first, second);
        std::cout << "distance " << aligned.distance << "\ncigar    " << aligned.cigar << "\n";
    } catch (const dinara::UnsupportedSymbols &error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
    return 0;
}
