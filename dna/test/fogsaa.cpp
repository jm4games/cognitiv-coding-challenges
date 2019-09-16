#include "fogsaa.hpp"
#include <limits>
#include <cstdint>
#include <queue>

using namespace std;

namespace dna
{

static const int64_t MatchScore = 1;
static const int64_t MisMatchScore = -1;
static const int64_t GapPenalty = -2;

static constexpr int64_t future_score_base(int64_t score, int64_t x1, int64_t x2)
{
    if (x2 < x1) {
        return x2 * score + GapPenalty * (x1 - x2);
    }

    return x1 * score + GapPenalty * (x2 - x1);
}

static constexpr int64_t fs_min(int64_t x1, int64_t x2)
{
    return future_score_base(MisMatchScore, x1, x2);
}

static constexpr int64_t fs_max(int64_t x1, int64_t x2)
{
    return future_score_base(MatchScore, x1, x2);
}

enum pairing_type : char
{
    Match,
    MisMatch,
    Gap
};

struct alignment
{
    int64_t ft_max = numeric_limits<int64_t>::min();
    pairing_type pairing = Match;
};

struct fitness_score
{
    int64_t min = 0;
    int64_t max = 0;
};

struct route
{
    fitness_score ft;
    int64_t s1_offset = 0;
    int64_t s2_offset = 0;

    bool operator() (const route& lhs, const route& rhs) const
    {
        if (lhs.ft.max < rhs.ft.max)
            return true;
        else if (lhs.ft.max == rhs.ft.max)
            return lhs.ft.min < rhs.ft.min;
        else
            return false;
    }
};

struct route_choice
{
    bool match;
    fitness_score non_gap;
    fitness_score gap_s1;
    fitness_score gap_s2;
};

static route_choice eval_route_choices(
        bool match, int64_t score, size_t s1_sz, size_t s2_sz, int64_t p1, int64_t p2)
{
    route_choice result;
    result.match = match;

    // compare
    int p1n = p1 + 1, p2n = p2 + 1;
    if (match) {
        int64_t score_nxt = score + MatchScore;
        int64_t x1 = s2_sz - p2n, x2 = s1_sz - p1n;
        result.non_gap.min = score_nxt + fs_min(x1, x2);
        result.non_gap.max = score_nxt + fs_max(x1, x2);
    } else {
        int64_t score_nxt = score + MisMatchScore;
        int64_t x1 = s2_sz - p2n, x2 = s1_sz - p1n;
        result.non_gap.min = score_nxt + fs_min(x1, x2);
        result.non_gap.max = score_nxt + fs_max(x1, x2);
    }

    // TODO: make sure gaps are corrrectly named

    // gap s1
    int64_t score_nxt = score + GapPenalty;
    int64_t x1 = s2_sz - p2n, x2 = s1_sz - p1;
    result.gap_s1.min = score_nxt + fs_min(x1, x2);
    result.gap_s1.max = score_nxt + fs_max(x1, x2);

    score_nxt = score + GapPenalty;
    x1 = s2_sz - p2, x2 = s1_sz - p1n;
    result.gap_s2.min = score_nxt + fs_min(x1, x2);
    result.gap_s2.max = score_nxt + fs_max(x1, x2);

    return result;
}

template<HelixStream T>
static constexpr void fill_helix_vector(T& helix, vector<byte> &vec)
{
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
fogsaa<T>::fogsaa(T& helix)
{
    fill_helix_vector(helix, s1_);
}

template<HelixStream T>
alignment_result fogsaa<T>::align_with(T& helix)
{
    vector<byte> s2;
    fill_helix_vector(helix, s2);

    // TODO: keep track of cur score, if t_max < cur score then bail

    // TODO: is it ok to assume no gaps in largest strand?
    vector<alignment> bestAlignment(max(s1_.size(), s2.size()), alignment{});
    priority_queue<route, vector<route>, route> priQueue;
    int64_t score = 0;
    int64_t p1 = 0, p2 = 0;

    //while (true) // has nodes
    //{
    //    int64_t x1 = s2.size() - p2, x2 = s1.size() - p1;
    //    int64_t fitMin = curScore + fs_min(x1, x2);
    //    int64_t fitMax = curScore + fs_max(x1, x2);
    if (p1 + 1 == s1_.size() || p2 + 1 == s2.size()) // leaf detection
    {
        // TODO: assign gap score to remainder of s1/s2.
    }

    route_choice choice = eval_route_choices(
            s1_[p1] == s2[p2], score, s1_.size(), s2.size(), p1, p2);
    if (choice.non_gap.max > choice.gap_s1.max && choice.non_gap.max > choice.gap_s2.max)
    {
    } else if (choice.gap_s1.max > choice.gap_s2.max)
    {
    }
    else
    {
    }
        // now pick one of comp, gap s1, or gap s2
    //}

    // TODO:
    return alignment_result{};
}

}
