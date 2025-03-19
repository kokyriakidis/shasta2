#pragma once

#include "cstdint.hpp"

namespace shasta {
    class Detangler;
    class TrivialDetangler;

    class Anchors;
    class AssemblyGraph;
    class Tangle;
}



class shasta::Detangler {
public:
    virtual bool operator()(const Anchors&, AssemblyGraph&, Tangle&) = 0;
};



class shasta::TrivialDetangler : public Detangler {
public:
    bool operator()(const Anchors&, AssemblyGraph&, Tangle&);
};
