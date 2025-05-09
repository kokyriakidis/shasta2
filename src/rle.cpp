#include "rle.hpp"
#include "Base.hpp"
using namespace shasta;

void shasta::rle(const vector<Base>& sequence, vector<Base>& rleSequence)
{
    rleSequence.clear();

    for(const Base& b: sequence) {
        if(rleSequence.empty() or b != rleSequence.back()) {
            rleSequence.push_back(b);
        }
    }
}
