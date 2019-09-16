#pragma once

#include "person.hpp"
#include "sequence_aligner.hpp"

namespace dna
{

template<HelixStream T>
class fogsaa {
    // TODO: explain why s1_ is here (caching!!)
    const std::vector<std::byte> s1_;
public:
    explicit fogsaa(T& helix);

    alignment_result align_with(T& helix);
};
}
