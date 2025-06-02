#pragma once

#include "Detangler.hpp"

#include "cstdint.hpp"

namespace shasta {
    class ExactDetangler;
}


// The ExactDetangler simply connects entrance/exit pairs that have a non-zero tangle matrix element.
// That is, they have one or more common reads.
// This has a tendency to generate breaks in the assembly.
class shasta::ExactDetangler : public Detangler {
public:
    bool operator()(Tangle&);
    bool operator()(Tangle3&);

};
