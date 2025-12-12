// Shasta.
#include "DisjointSets.hpp"
#include "orderPairs.hpp"
using namespace shasta2;

// Standard library.
#include "iostream.hpp" // Remove when done debugging.
#include "utility.hpp"


DisjointSets::DisjointSets(uint64_t n) :
    rank(n),
    parent(n),
    disjointSets(&rank[0], &parent[0])
{
    initializeDisconnected();
}



void DisjointSets::initializeDisconnected()
{
    const uint64_t n = rank.size();
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }
}



void DisjointSets::initializeDisconnected(uint64_t i)
{
    disjointSets.make_set(i);
}



void DisjointSets::unionSet(uint64_t i, uint64_t j)
{
    disjointSets.union_set(i, j);
}



void DisjointSets::link(uint64_t i, uint64_t j)
{
    disjointSets.link(i, j);
}



uint64_t DisjointSets::findSet(uint64_t i)
{
    return disjointSets.find_set(i);
}


// Find all components of size at least equal to minComponentSize
// and gather them sorted by decreasing size.
void DisjointSets::gatherComponents(
    uint64_t minComponentSize,
    vector< vector<uint64_t> >& components
    )
{
    // Gather all components, including the dummy empty ones.
    const uint64_t n = rank.size();
    vector< vector<uint64_t> > allComponents(n);
    for(uint64_t i=0; i<n; i++) {
        const uint64_t componentId = findSet(i);
        allComponents[componentId].push_back(i);
    }

    // Create a table with the indexes of the ones larger than minComponentSize
    // and their size.
    vector< pair<uint64_t, uint64_t> > componentTable;
    for(uint64_t i=0; i<n; i++) {
        const vector<uint64_t>& component = allComponents[i];
        if(component.size() >= minComponentSize) {
            componentTable.emplace_back(i, component.size());
        }
    }

    // Sort the componentTable by decreasing component size.
    std::ranges::sort(componentTable, OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    // Store the components in this order.
    components.clear();
    for(const auto& p: componentTable) {
        const uint64_t i = p.first;
        const vector<uint64_t>& component = allComponents[i];
        components.emplace_back(component);
    }
}

