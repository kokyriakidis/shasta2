#pragma once

#include <iosfwd.hpp>
#include "vector.hpp"

namespace shasta {

    class AlignedBase;
    void simpleFastaRead(istream&, vector< vector<AlignedBase> >&);
}
