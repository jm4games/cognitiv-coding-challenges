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
            return;
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
    int64_t p1 = 0, p2 = 0;
    vector<byte> s2;
    vector<byte> bestAlignment;
    vector<byte> curAlignment;
    //queue<..> priQueue; ordered by fitMax
    fillHelixVector(helix, s2);

    int64_t curScore = 0;

    //while (true) // has nodes
    //{
    //    int64_t x1 = s2.size() - p2, x2 = s1.size() - p1;
    //    int64_t fitMin = curScore + fsMin(x1, x2);
    //    int64_t fitMax = curScore + fsMax(x1, x2);

    //    // compare
    //    int p1n = p1 + 1, p2n = p2 + 1;
    //    if (s1[p1n] == s2[p2n]) {
    //        int64_t scoreNxt = curScore + MatchScore;
    //        int64_t x1 = s2.size() - p2n, x2 = s1.size() - p1n;
    //        int64_t fitMinN = scoreNxt + fsMin(x1, x2);
    //        int64_t fitMaxN = scoreNxt + fsMax(x1, x2);
    //    } else {
    //        int64_t scoreNxt = curScore + MisMatchScore;
    //        int64_t x1 = s2.size() - p2n, x2 = s1.size() - p1n;
    //        int64_t fitMinN = scoreNxt + fsMin(x1, x2);
    //        int64_t fitMaxN = scoreNxt + fsMax(x1, x2);
    //    }

    //    // gap s1
    //    int64_t scoreNxt = curScore + GapPenalty;
    //    int64_t x1 = s2.size() - p2n, x2 = s1.size() - p1;
    //    int64_t fitMinN = scoreNxt + fsMin(x1, x2);
    //    int64_t fitMaxN = scoreNxt + fsMax(x1, x2);

    //    // gap s2
    //    int64_t scoreNxt = curScore + GapPenalty;
    //    int64_t x1 = s2.size() - p2, x2 = s1.size() - p1n;
    //    int64_t fitMinN = scoreNxt + fsMin(x1, x2);
    //    int64_t fitMaxN = scoreNxt + fsMax(x1, x2);

        // now pick one of comp, gap s1, or gap s2
    //}

    // TODO:
    return AlignmentReport{};
}

}
