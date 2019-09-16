#include "person_analyzer.hpp"

namespace dna
{

analyzer::analyzer() { }

template <typename P> requires (Person<P> && IsPerson<P>)
     int analyzer::analyze_people_async(P& p1, P& p2)
{
    return 0;
}

} // dna
