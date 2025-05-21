#pragma once

namespace shasta {
    class Detangler;
    class Tangle;
    class Tangle2;
    class Tangle3;
}

class shasta::Detangler {
public:
    virtual bool operator()(Tangle&) = 0;
    virtual bool operator()(Tangle2&) = 0;
    virtual bool operator()(Tangle3&) = 0;

    bool debug = false;
};
