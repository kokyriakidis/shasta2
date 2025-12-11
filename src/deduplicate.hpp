#pragma once

// Shasta.
#include "algorithm.hpp"
#include "SHASTA2_ASSERT.hpp"

// Standard library.
#include "iostream.hpp"
#include "span.hpp"
#include "vector.hpp"

namespace shasta2 {

    // Remove duplicate elements in a vector.
    template<class T> void deduplicate(vector<T>& v)
    {
        sort(v.begin(), v.end());
        v.resize(unique(v.begin(), v.end()) - v.begin());
    }


    // Remove duplicate elements in a vector and count occurrences of each
    template<class T, class Int> void deduplicateAndCount(
        vector<T>& v,
        vector<Int>& count,
        bool alreadySorted = false)
    {
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
        typename vector<T>::iterator output = v.begin();
        typename vector<T>::iterator input = v.begin();
        while(input != v.end()) {

            // Store this element.
            *output = *input;
            ++output;

            // Count how many there are.
            typename vector<T>::iterator it = input;
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



    // Remove duplicate elements in a vector and count occurrences of each.
    // Keep only the ones that occur at least minCount times.
    template<class T, class Int> void deduplicateAndCountWithThreshold(
        vector<T>& v,
        vector<Int>& count,
        Int minCount
        )
    {
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
        typename vector<T>::iterator output = v.begin();
        typename vector<T>::iterator input = v.begin();
        while(input != v.end()) {


            // Count how many there are.
            typename vector<T>::iterator it = input;
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



    // Remove duplicate elements in a span and count occurrences of each.
    // Keep only the ones that occur at least minCount times.
    template<class T, class Int> void deduplicateAndCountWithThreshold(
        span<T>& v,
        vector<Int>& count,
        Int minCount
        )
    {
        // Clear the count vector.
        count.clear();

        // If the given vector is empty, return now.
        if(v.empty()) {
            return;
        }

        // Sort the span.
        sort(v.begin(), v.end());

        // Add elements, keeping track of the number
        // of occurrences of each.
        typename span<T>::iterator output = v.begin();
        typename span<T>::iterator input = v.begin();
        while(input != v.end()) {


            // Count how many there are.
            typename span<T>::iterator it = input;
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
    }



    // Remove duplicate elements in a vector and count occurrences of each.
    // Keep only the ones that occur exactly once.
    template<class T> void deduplicateAndCountAndKeepUnique(
        vector<T>& v)
    {

        // If the given vector is empty, return now.
        if(v.empty()) {
            return;
        }

        // Sort the vector.
        sort(v.begin(), v.end());

        // Add elements, keeping track of the number
        // of occurrences of each.
        typename vector<T>::iterator output = v.begin();
        typename vector<T>::iterator input = v.begin();
        while(input != v.end()) {


            // Count how many there are.
            typename vector<T>::iterator it = input;
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



    // Remove duplicate elements in a vector and count occurrences of each.
    // Keep only the ones that occur exactly once.
    // Version that uses a custom comparator.
    template<class T, class Comparator> void deduplicateAndCountAndKeepUnique(
        vector<T>& v,
        const Comparator& comparator)
    {

        // If the given vector is empty, return now.
        if(v.empty()) {
            return;
        }

        // Sort the vector.
        sort(v.begin(), v.end(), comparator);

        // Add elements, keeping track of the number
        // of occurrences of each.
        typename vector<T>::iterator output = v.begin();
        typename vector<T>::iterator input = v.begin();
        while(input != v.end()) {


            // Count how many there are.
            typename vector<T>::iterator it = input;
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

