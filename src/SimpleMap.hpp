#pragma once

/******************************************************************

A simple map with integer keys.

It uses open addressing and linear probing.
See https://en.wikipedia.org/wiki/Open_addressing

This uses the lowest order bits of the keys as a hash.
Therefore it is suited if the keys are "random".

The Value type must be default constructible.

******************************************************************/

#include "invalid.hpp"
#include "SHASTA2_ASSERT.hpp"

#include "cstdint.hpp"
#include <type_traits>
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {
    template <class Key, class Value> class SimpleMap;
    void testSimpleMap();
}


template <class Key, class Value> class shasta::SimpleMap
{
public:

    // The Key is required to be an unsigned integer type.
    static_assert(std::is_integral<Key>());
    static_assert(std::is_unsigned<Key>());

    using value_type = pair<Key, Value>;

    // A vector that will hold (Key, Value) pairs.
    // Empty slots have the Key set to invalid<Key>.
    vector<value_type> data;
    Key mask;
    uint64_t itemCount = 0UL;

    // Construct with n slots.
    // The actual number of slots is rounded up to the next power of 2.
    SimpleMap(uint64_t n) :
        data(1UL << (64 - __builtin_clzl(n)), value_type(invalid<Key>, Value())),
        mask(data.size() - 1UL)
    {}



    // Adds a new (Key, value) pair.
    // Returns a pointer to the newly added pair.
    // If the key already exists, this asserts.
    value_type* insertNew(const value_type& p)
    {
        SHASTA2_ASSERT(itemCount < data.size());

        const Key& key = p.first;
        SHASTA2_ASSERT(key != invalid<Key>);

        uint64_t slotIndex = key;
        while(true) {
            slotIndex &= mask;
            value_type& slot = data[slotIndex];
            const Key& slotKey = slot.first;
            SHASTA2_ASSERT(slotKey != key);
            if(slot.first == invalid<Key>) {
                slot = p;
                ++itemCount;
                return &slot;
            }
            ++slotIndex;
        }
    }



    // Returns a pointer to an existing (Key, Value) pair.
    // If the key does not exists, this returns 0.
    value_type* getExisting(const Key& key)
    {
        SHASTA2_ASSERT(itemCount < data.size());
        SHASTA2_ASSERT(key != invalid<Key>);

        uint64_t slotIndex = key;
        while(true) {
            slotIndex &= mask;
            value_type& slot = data[slotIndex];
            const Key& slotKey = slot.first;
            if(slotKey == invalid<Key>) {
                return 0;
            }
            if(slotKey == key) {
                return &slot;
            }
            ++slotIndex;
        }
    }


    // If the key exists, this returns a pointer to the (Key, Value)
    // pair with that Key and does not modify the Value
    // If the key does nto exist, this adds a new (Key, Value) pair.
    // and returns a pointer to the newly added pair.
    value_type* insertNewOrGetExisting(const value_type& p)
    {
        SHASTA2_ASSERT(itemCount < data.size());

        const Key& key = p.first;
        SHASTA2_ASSERT(key != invalid<Key>);

        uint64_t slotIndex = key;
        while(true) {
            slotIndex &= mask;
            value_type& slot = data[slotIndex];
            const Key& slotKey = slot.first;
            if(slotKey == key) {
                return &slot;
            }
            if(slot.first == invalid<Key>) {
                slot = p;
                ++itemCount;
                return &slot;
            }
            ++slotIndex;
        }
    }
};
