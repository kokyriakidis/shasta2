#include "SimpleMap.hpp"
using namespace shasta2;

#include "iostream.hpp"



void shasta2::testSimpleMap()
{
    class A {
    public:
        double x = 3.5;
    };

    const uint64_t n = 10;
    SimpleMap<uint64_t, A> m(n);
    cout << "Mask is " << m.mask << endl;

    m.insertNew(make_pair(102, A(102.)));
    m.insertNew(make_pair(118, A(118.)));
    m.insertNew(make_pair(10, A(10.)));

    for(uint64_t i=0; i<m.data.size(); i++) {
        const auto& p = m.data[i];
        cout << i << " " << p.first << " " << p.second.x << endl;
    }

    const auto p = m.getExisting(102);
    if(p) {
        cout << p->first << " " << p->second.x << endl;
    } else {
        cout << "Not found." << endl;
    }
}
