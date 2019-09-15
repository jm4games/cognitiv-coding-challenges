#pragma once

#include "person.hpp"

namespace dna
{

struct AlignmentReport {
};

template<HelixStream T>
class Fogsaa {

public:
    AlignmentReport runGlobalSequenceAlignment(T& a, T& b);
};
}
