#pragma once

#include "sequence_buffer.hpp"
#include "person.hpp"
#include <cstdint>
#include <vector>

namespace dna {

// Defines the location of a subsequence in a helix.
struct location
{
    int64_t offset = 0;
    int64_t size = 0;
};

// Defines the location where two strands are conflicting.
struct mutation
{
    location helix1;
    location helix2;
};

// Describes how similar two helixes are and where they differ (mutations).
// If aligment fails (eg due to some side effect) the error value will be set
// to a non-empty value.
struct alignment_result
{
    std::vector<mutation> mutations;
    double simularityScore = 0;
    std::string error;
};

template<HelixStream T>
class SequenceAligner {
public:
    virtual alignment_result align_with(T& a, T& b) = 0;
};

} // dna
