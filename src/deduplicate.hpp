#pragma once

// Shasta.
#include "algorithm.hpp"
#include "SHASTA2_ASSERT.hpp"

// Standard library.
#include "iostream.hpp"
#include "span.hpp"
#include "vector.hpp"

namespace shasta2 {

    // Remove duplicate elements in a vector-like container.
    template<class Vector> void deduplicate(Vector& v)
    {
        sort(v.begin(), v.end());
        v.resize(unique(v.begin(), v.end()) - v.begin());
    }


    // Remove duplicate elements in a vector-like and count occurrences of each
    template<class Vector, class Int> void deduplicateAndCount(
        Vector& v,
        vector<Int>& count,
        bool alreadySorted = false)
    {
        using iterator = typename Vector::iterator;

        // Clear the count vector.
        count.clear();

        // If the given vector is empty, return now.
        if(v.empty()) {
            return;
        }

        // Sort the vector.
        if(not alreadySorted) {
            sort(v.begin(), v.end());
        }

        // Add elements, keeping track of the number
        // of occurrences of each.
        iterator output = v.begin();
        iterator input = v.begin();
        while(input != v.end()) {

            // Store this element.
            *output = *input;
            ++output;

            // Count how many there are.
            iterator it = input;
            while(it!=v.end() && *it==*input) {
                ++it;
            }

            // Store the count.
            count.push_back(Int(it - input));

            // Update our output iterator.
            input = it;

        }
        v.resize(count.size());
    }



    // Remove duplicate elements in a vector-like and count occurrences of each.
    // Keep only the ones that occur at least minCount times.
    template<class Vector, class Int> void deduplicateAndCountWithThreshold(
        Vector& v,
        vector<Int>& count,
        Int minCount
        )
    {
        using iterator = typename Vector::iterator;

        // Clear the count vector.
        count.clear();

        // If the given vector is empty, return now.
        if(v.empty()) {
            return;
        }

        // Sort the vector.
        sort(v.begin(), v.end());

        // Add elements, keeping track of the number
        // of occurrences of each.
        iterator output = v.begin();
        iterator input = v.begin();
        while(input != v.end()) {


            // Count how many there are.
            iterator it = input;
            while(it!=v.end() && *it==*input) {
                ++it;
            }
            const Int n = Int(it - input);

            if(n >= minCount) {

                // Store this element.
                *output = *input;
                ++output;

                // Store the count.
                count.push_back(n);
            }

            // Update our output iterator.
            input = it;

        }
        v.resize(count.size());
    }



    // Remove duplicate elements in a vector-like and count occurrences of each.
    // Keep only the ones that occur exactly once.
    template<class Vector> void deduplicateAndCountAndKeepUnique(
        Vector& v)
    {
        using iterator = typename Vector::iterator;

        // If the given vector is empty, return now.
        if(v.empty()) {
            return;
        }

        // Sort the vector.
        sort(v.begin(), v.end());

        // Add elements, keeping track of the number
        // of occurrences of each.
        iterator output = v.begin();
        iterator input = v.begin();
        while(input != v.end()) {


            // Count how many there are.
            iterator it = input;
            while(it!=v.end() && *it==*input) {
                ++it;
            }
            const uint64_t n = it - input;

            if(n == 1) {

                // Store this element.
                *output = *input;
                ++output;
            }

            // Update our output iterator.
            input = it;

        }
        v.resize(output - v.begin());
    }



    // Remove duplicate elements in a vector-like and count occurrences of each.
    // Keep only the ones that occur exactly once.
    // Version that uses a custom comparator.
    template<class Vector, class Comparator> void deduplicateAndCountAndKeepUnique(
        Vector& v,
        const Comparator& comparator)
    {
        using iterator = typename Vector::iterator;

        // If the given vector is empty, return now.
        if(v.empty()) {
            return;
        }

        // Sort the vector.
        sort(v.begin(), v.end(), comparator);

        // Add elements, keeping track of the number
        // of occurrences of each.
        iterator output = v.begin();
        iterator input = v.begin();
        while(input != v.end()) {


            // Count how many there are.
            iterator it = input;
            while(it!=v.end() and (not comparator(*it, *input)) and (not comparator(*input, *it))) {
                ++it;
            }
            const uint64_t n = it - input;

            if(n == 1) {

                // Store this element.
                *output = *input;
                ++output;
            }

            // Update our output iterator.
            input = it;

        }
        v.resize(output - v.begin());
    }



    inline void testDeduplicateAndCount()
    {
        vector<int> v = {7, 4, 5, 7, 4, 18, 2, 4};
        vector<int> count;
        deduplicateAndCountWithThreshold(v, count, 2);
        SHASTA2_ASSERT(v.size() == count.size());
        for(uint64_t i=0; i<v.size(); i++) {
            cout << v[i] << " " << count[i] << endl;
        }
    }
}

