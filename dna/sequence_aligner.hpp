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
    int64_t length = 0;

    bool operator==(const location& l) const
    {
        return offset == l.offset && length == l.length;
    }
};

// Defines the location where two strands are conflicting.
struct mutation
{
    location helix1;
    location helix2;

    mutation(location h1, location h2): helix1(h1), helix2(h2) {}

    bool operator==(const mutation& m) const
    {
        return helix1 == m.helix1 && helix2 == m.helix2;
    }
};

// Describes how similar two helixes are and where they differ (mutations).
// If aligment fails (eg due to some side effect) the error value will be set
// to a non-empty value.
struct alignment_result
{
    std::vector<mutation> mutations;
    double simularity_score = 0;
    std::string error;
};

template<HelixStream T>
class sequence_aligner {
public:
    virtual alignment_result align(T& a, T& b) const = 0;
};

} // dna
