#pragma once


// BitCounter is a traits class that can be used to define the
// number of bits in an integer type.
// std::numeric_limits::digits only works for built-in types
// and cannot be extended.

// Standard library.
#include "cstdint.hpp"

namespace shasta {

    template<class Int> class BitCounter {
    public:
        static constexpr int numberOfBits = 0;
    };

    template<> class BitCounter<uint8_t> {
    public:
        static constexpr int numberOfBits = 8;
    };

    template<> class BitCounter<uint16_t> {
    public:
        static constexpr int numberOfBits = 16;
    };

    template<> class BitCounter<uint32_t> {
    public:
        static constexpr int numberOfBits = 32;
    };

    template<> class BitCounter<uint64_t> {
    public:
        static constexpr int numberOfBits = 64;
        using doubleSizeType = __uint128_t;
    };

    template<> class BitCounter<__uint128_t> {
    public:
        static constexpr int numberOfBits = 128;
    };
}
