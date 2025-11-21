#pragma once


#include <Base.hpp>

#include "vector.hpp"

// Convert a vector<Base> to RLE.
namespace shasta {
    void rle(const vector<Base>& sequence, vector<Base>& rleSequence);
}
