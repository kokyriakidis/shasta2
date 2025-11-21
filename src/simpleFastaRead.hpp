#pragma once

#include "iosfwd.hpp"
#include "vector.hpp"

namespace shasta2 {

    class AlignedBase;
    void simpleFastaRead(istream&, vector< vector<AlignedBase> >&);
}
