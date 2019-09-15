#include "fogsaa.hpp"
#include <cstdint>

using namespace std;

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
static constexpr void fillHelixVector(T& helix, vector<byte> &vec) {
    vec.reserve(helix.size());
    while (true) {
        auto buf = helix.read();
        if (buf.size() == 0) {
            break;
        }

        auto end = buf.end();
        for (auto it = buf.begin(); it != end; ++it) {
            vec.emplace_back(*it);
        }
    }
}

template<HelixStream T>
Fogsaa<T>::Fogsaa(T& helix)
{
    fillHelixVector(helix, s1);
}

template<HelixStream T>
AlignmentReport Fogsaa<T>::alignWith(T& helix)
{
    vector<byte> s2;
    fillHelixVector(helix, s2);

    // TODO:
    return AlignmentReport{};
}

}
