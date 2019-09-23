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
    std::string error;
    double similarity_score = 0;

    alignment_result() {}

    alignment_result(std::vector<mutation>&& muts, std::string&& err, double score)
        : mutations(std::move(muts)), error(std::move(err)), similarity_score(score)
    {}

    alignment_result(const alignment_result& other)
        : mutations(other.mutations),
          similarity_score(other.similarity_score)
    {
        if (other.error != "")
            error = std::move(other.error);
    }

    alignment_result(alignment_result&& other) noexcept
        : mutations(std::move(other.mutations)),
          error(std::move(other.error)),
          similarity_score(other.similarity_score)
    {}

    alignment_result& operator=(const alignment_result& other)
    {
        mutations = other.mutations;
        error = other.error;
        similarity_score = other.similarity_score;
        return *this;
    }

    alignment_result& operator=(alignment_result&& other)
    {
        mutations = std::move(other.mutations);
        if (other.error != "")
            error = std::move(other.error);
        similarity_score = other.similarity_score;
        return *this;
    }
};

template<HelixStream T>
class sequence_aligner {
public:
    virtual alignment_result align(T& a, T& b) const = 0;
};

} // dna
