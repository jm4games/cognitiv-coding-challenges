#include "fogsaa.hpp"
#include <cstdint>

namespace dna
{

static const int64_t MatchScore = 1;
static const int64_t MisMatchScore = -1;
static const int64_t GapPenalty = -2;

static constexpr int64_t futureScoreBase(int64_t score, int64_t x1, int64_t x2)
{
    if (x2 < x1) {
        return x2 * score + GapPenalty * (x1 - x2);
    }

    return x1 * score + GapPenalty * (x2 - x1);
}

static constexpr int64_t fsMin(int64_t x1, int64_t x2)
{
    return futureScoreBase(MisMatchScore, x1, x2);
}

static constexpr int64_t fsMax(int64_t x1, int64_t x2)
{
    return futureScoreBase(MatchScore, x1, x2);
}

template<HelixStream T>
AlignmentReport Fogsaa<T>::runGlobalSequenceAlignment(T& a, T& b)
{
    // TODO:
    return AlignmentReport{};
}

}
