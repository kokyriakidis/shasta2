#pragma once


namespace shasta2 {
    template<class RandomIterator> inline bool isPeriodic(RandomIterator begin, RandomIterator end, uint64_t period);
}



// Return true if the sequence defined by the given iterators is periodic with period p.
template<class RandomIterator> inline bool shasta2::isPeriodic(
    RandomIterator begin, RandomIterator end, uint64_t period)
{
    if((end - begin) % period) {
        return false;
    }

    for(RandomIterator itA=begin; ; ++itA) {
        const RandomIterator itB = itA + period;
        if(itB >= end) {
            break;
        }
        if(*itA != *itB) {
            return false;
        }
    }
    return true;
}
