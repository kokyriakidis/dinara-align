/* This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
 * MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/. */
// dinara-align from C++: the edit distance and an optimal alignment's CIGAR for two sequences, and the
// query placed inside the reference under gap-affine costs.
//
//     pixi run build-c
//     c++ -std=c++17 -I build/c c/example.cpp -L build/c -ldinara -Wl,-rpath,build/c -o example
//     ./example ACGTACGTTTGCA ACGTCGTTTTGCA
#include <iostream>

#include "dinara.h"

// Aligns the two sequences given, or two built in, globally under unit costs and then the query inside the
// reference under gap-affine costs.
int main(int argc, char **argv) {
    std::string_view reference = argc > 2 ? argv[1] : "ACGTACGTTTGCA";
    std::string_view query = argc > 2 ? argv[2] : "ACGTCGTTTTGCA";
    try {
        dinara::Alignment aligned = dinara::align(reference, query);
        std::cout << "distance " << aligned.cost << "\ncigar    " << aligned.cigar << "\n";
        dinara::Alignment placed = dinara::align(reference, query, dinara::Costs::affine(4, 6, 2), dinara::Mode::infix());
        std::cout << "placed   " << placed.cigar << " at " << placed.reference_start << ".." << placed.reference_end
                  << ", cost " << placed.cost << "\n";
    } catch (const std::invalid_argument &error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
    return 0;
}
