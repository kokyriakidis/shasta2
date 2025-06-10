#pragma once

namespace shasta {
    class Detangler;
    class Tangle;
}

class shasta::Detangler {
public:

    // This evaluates the TangleMatrix and:
    // - If detangling is not possible, returns false and stores nothing in the Tangle.
    // - If detangling is possible:
    //   * It calls Tangle::connect to store in the Tangle the entrance/exit pairs to be connected.
    //   * If doDetangle is true, it calls Tangle::detangle.
    //   * It returns true.
    virtual bool operator()(Tangle&, bool doDetangle) = 0;

    bool debug = false;
};
