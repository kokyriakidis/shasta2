#pragma once

namespace shasta2 {
    class Detangler;
    class Tangle1;
}

class shasta2::Detangler {
public:

    // This evaluates the TangleMatrix and:
    // - If detangling is not possible, returns false and stores nothing in the Tangle.
    // - If detangling is possible:
    //   * It calls Tangle::connect to store in the Tangle the entrance/exit pairs to be connected.
    //   * It calls Tangle::detangle.
    //   * It returns true.
    virtual bool operator()(Tangle1&)
    {
        return false;
    }

    bool debug = false;
};
