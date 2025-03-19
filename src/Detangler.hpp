#pragma once

#include "cstdint.hpp"

namespace shasta {
    class Detangler;
    class TrivialDetangler;

    class Tangle;
}



class shasta::Detangler {
public:
    virtual bool operator()(Tangle&) = 0;
};



class shasta::TrivialDetangler : public Detangler {
public:
    bool operator()(Tangle&);
};
