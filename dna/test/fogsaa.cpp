#include "fogsaa.hpp"
#include <limits>
#include <cstdint>
#include <queue>
#include <algorithm>

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
    None,
    Match,
    MisMatch,
    Gap
};

static constexpr int64_t pairing_score(pairing_type type)
{
    switch(type)
    {
        case Match:
            return MatchScore;
        case MisMatch:
            return MisMatchScore;
        case Gap:
            return GapPenalty;
        default:
            return 0;
    }
}

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

struct pairing
{
    fitness_score ft;
    int64_t s1_offset = -1;
    int64_t s2_offset = -1;
    pairing_type type = None;

    bool operator() (const pairing& lhs, const pairing& rhs) const
    {
        if (lhs.ft.max < rhs.ft.max)
            return true;
        else if (lhs.ft.max == rhs.ft.max)
            return lhs.ft.min < rhs.ft.min;
        else
            return false;
    }

    bool is_set() const
    {
        return type == None;
    }
};

struct pairing_choice
{
    bool match = false;
    pairing non_gap;
    pairing gap_s1;
    pairing gap_s2;
};

struct pairing_key
{
    int64_t p1 = 0;
    int64_t p2 = 0;

    bool operator==(const pairing_key other) const
    {
        return p1 == other.p1 && p2 == other.p2;
    }
};

struct hash_pairing_key
{
    size_t operator() (const pairing_key& k) const {
        return hash<int64_t>()(k.p1) ^ (hash<int64_t>()(k.p2) << 1);
    }
};

static pairing_choice eval_pairing_choices(
        const vector<byte>& s1, const vector<byte>& s2, int64_t score, int64_t p1, int64_t p2)
{
    int64_t s1_sz = s1.size(), s2_sz = s2.size();
    pairing_choice result;

    // compare
    int64_t p1n = p1 + 1, p2n = p2 + 1;
    result.match = s1[p1n] == s2[p2n];

    int64_t x1 = s2_sz - p2n, x2 = s1_sz - p1n;
    int64_t score_nxt = score + (result.match ? MatchScore : MisMatchScore);
    result.non_gap.s1_offset = p1n;
    result.non_gap.s2_offset = p2n;
    result.non_gap.ft = fitness_score{score_nxt + fs_min(x1, x2),  score_nxt + fs_max(x1, x2)};

    score_nxt = score + GapPenalty;

    // gap s1
    x1 = s2_sz - p2n, x2 = s1_sz - min(p1, (int64_t)0);
    result.gap_s1.s1_offset = p1;
    result.gap_s1.s2_offset = p2n;
    result.gap_s1.ft = fitness_score{score_nxt + fs_min(x1, x2), score_nxt + fs_max(x1, x2)};

    // gap s2
    x1 = s2_sz - min(p2, (int64_t)0), x2 = s1_sz - p1n;
    result.gap_s2.s1_offset = p1n;
    result.gap_s2.s2_offset = p2;
    result.gap_s2.ft = fitness_score{score_nxt + fs_min(x1, x2), score_nxt + fs_max(x1, x2)};

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

    vector<alignment> bestAlignment(max(s1_.size(), s2.size()), alignment{});
    vector<alignment> curAlignment(max(s1_.size(), s2.size()), alignment{});
    unordered_map<pairing_key, int64_t, hash_pairing_key> knownPairings;
    priority_queue<pairing, vector<pairing>, pairing> priQueue;
    int64_t score = 0, best_score = 0, p1 = -1, p2 = -1;
    pairing cur_pairing;

    pairing_choice choice = eval_pairing_choices(s1_, s2, score, p1, p2);
    priQueue.push(choice.non_gap);
    priQueue.push(choice.gap_s1);
    priQueue.push(choice.gap_s2);

    // TODO: keep track of cur score, if t_max < cur score then bail

    while (!priQueue.empty())
    {
        if (!cur_pairing.is_set()) {
            cur_pairing = priQueue.pop();

            if (cur_pairing.ft.max <= best_score) {
                cur_pairing = pairing{}; // reset so we can pop
                continue;
            }
        }

        score += pairing_score(cur_pairing.type);
        if (cur_pairing.s1_offset + 1 == s1_.size() || cur_pairing.s2_offset + 1 == s2.size()) // leaf detection
        {
            // TODO: assign gap score to remainder of s1/s2.
            // TODO: update best pairing if cur pairing better
            continue;
        }

        // TODO: commit current pairing once we have child pairing that is better then existing
        // TODO: compare current best pairing with existing best pairing, if its worse we need to back
        // track
        pairing_choice choice = eval_pairing_choices(s1_, s2, score, p1, p2);
        if (choice.non_gap.ft.max > choice.gap_s1.ft.max && choice.non_gap.ft.max > choice.gap_s2.ft.max)
        {
            auto existing = knownPairings.find(
                    pairing_key{choice.non_gap.s1_offset, choice.non_gap.s2_offset});
            if (existing != knownPairings.end()) {
                if (choice.non_gap.ft.max <= existing->second) {
                    // TODO: backtrack
                    continue;
                }
            }

            curAlignment.emplace_back(choice.non_gap.ft.max, Match);
            cur_pairing = choice.non_gap;

            priQueue.push(choice.gap_s1);
            priQueue.push(choice.gap_s2);

        } else if (choice.gap_s1.ft.max > choice.gap_s2.ft.max)
        {
        }
        else
        {
        }
    }

    return alignment_result{};
}

} // dna
