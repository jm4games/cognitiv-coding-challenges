#pragma once

#include "person.hpp"
#include "sequence_aligner.hpp"

#include <future>

namespace dna {

template<typename T>
concept bool ChromesomeSimliar = requires(T a, T b) {
    a.chromosomes() == b.chromosomes();
};

template<typename T>
concept bool IsPerson = requires(T a) {
    a.chromosomes() == 24;
};

class analyzer {

public:
    analyzer() {};

    template <typename P> requires (Person<P> && IsPerson<P>)
         //std::future<std::array<AlignmentReport, 23>> analyze_people_async(P& p1, P& p2);
         int analyze_people_async(P& p1, P& p2) { return 0; }
};

} // dna
