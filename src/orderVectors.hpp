#pragma once

#include "vector.hpp"

// Classes to sort vectors by size.

namespace shasta2 {

    template<class T> class OrderVectorsByIncreasingSize;
    template<class T> class OrderVectorsByDecreasingSize;

}



template<class T> class shasta2::OrderVectorsByIncreasingSize {
public:
     bool operator()(const vector<T>& x, const vector<T>& y) const
    {
         return x.size() < y.size();
    }
};



template<class T> class shasta2::OrderVectorsByDecreasingSize {
public:
     bool operator()(const vector<T>& x, const vector<T>& y) const
    {
         return x.size() > y.size();
    }
};

