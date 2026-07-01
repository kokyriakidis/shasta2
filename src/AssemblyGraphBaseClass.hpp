#pragma once

#include <boost/graph/adjacency_list.hpp>


namespace shasta2 {

    class AssemblyGraph;
    class AssemblyGraphVertex;
    class AssemblyGraphEdge;

    using AssemblyGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AssemblyGraphVertex,
        AssemblyGraphEdge>;
}
