#include "weightedShuffle.hpp"
using namespace shasta2;

#include "iostream.hpp"
#include "iterator.hpp"


void shasta2::testWeightedShuffle()
{

    const vector<double> weights = {1., 100., 3., 50., 10.};
    vector<uint64_t> shuffle;
    std::mt19937 generator;

    for(uint64_t iteration=0; iteration<30; iteration++) {
        weightedShuffle(weights, generator, shuffle);
        std::ranges::copy(shuffle, ostream_iterator<uint64_t>(cout, " "));
        cout << endl;
    }
}
