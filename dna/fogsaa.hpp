#pragma once

#include "person.hpp"
#include <vector>

namespace dna
{

struct AlignmentReport {
};

template<HelixStream T>
class Fogsaa {
    const std::vector<std::byte> s1;
public:
    explicit Fogsaa(T& helix);

    AlignmentReport alignWith(T& helix);
};
}
