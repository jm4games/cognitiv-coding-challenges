#pragma once

#include "person.hpp"
#include "sequence_aligner.hpp"

namespace dna
{

template<HelixStream T>
class Fogsaa {
    const std::vector<std::byte> s1_;
public:
    explicit Fogsaa(T& helix);

    AlignmentReport alignWith(T& helix);
};
}
