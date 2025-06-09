#pragma once

namespace shasta {
    class Detangler;
    class Tangle;
}

class shasta::Detangler {
public:
    virtual bool operator()(Tangle&) = 0;

    bool debug = false;
};
