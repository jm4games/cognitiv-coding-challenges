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
    Gap,
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

struct final_pairing
{
    byte s1 = static_cast<byte>(0);
    byte s2 = static_cast<byte>(0);

    final_pairing(byte b1, byte b2) : s1(b1), s2(b2) {}
};

struct fitness_score
{
    int64_t min = 0;
    int64_t max = 0;
};

struct pairing
{
    fitness_score ft;
    int64_t s1offset = -1;
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
    int64_t s1sz = s1.size(), s2_sz = s2.size();
    pairing_choice result;

    // compare
    int64_t p1n = p1 + 1, p2n = p2 + 1;
    result.non_gap.type = s1[p1n] == s2[p2n] ? Match : MisMatch;

    int64_t x1 = s2_sz - p2n, x2 = s1sz - p1n;
    int64_t score_nxt = score + (result.non_gap.type == Match ? MatchScore : MisMatchScore);
    result.non_gap.s1offset = p1n;
    result.non_gap.s2_offset = p2n;
    result.non_gap.ft = fitness_score{score_nxt + fs_min(x1, x2),  score_nxt + fs_max(x1, x2)};

    score_nxt = score + GapPenalty;

    // gap s1
    //x1 = s2_sz - p2n, x2 = s1sz - min(p1, (int64_t)0);
    //result.gap_s1.s1offset = p1;
    //result.gap_s1.s2_offset = p2n;
    //result.gap_s1.ft = fitness_score{score_nxt + fs_min(x1, x2), score_nxt + fs_max(x1, x2)};
    //result.gap_s1.type = GapS1;

    // gap s2
    x1 = s2_sz - min(p2, (int64_t)0), x2 = s1sz - p1n;
    result.gap_s2.s1offset = p1n;
    result.gap_s2.s2_offset = p2;
    result.gap_s2.ft = fitness_score{score_nxt + fs_min(x1, x2), score_nxt + fs_max(x1, x2)};
    result.gap_s1.type = Gap;

    return result;
}

alignment_result fogsaa::align_bytes(const vector<byte>& s1, const vector<byte>& s2) const
{
    vector<final_pairing> best_pairings;
    vector<final_pairing> cur_pairings;
    unordered_map<pairing_key, int64_t, hash_pairing_key> knownPairings;
    priority_queue<pairing, vector<pairing>, pairing> priQueue;
    int64_t score = 0, best_score = 0, base_offset = 0, p1 = -1, p2 = -1;
    pairing cur_pairing;

    best_pairings.reserve(max(s1.size(), s2.size()));
    cur_pairings.reserve(max(s1.size(), s2.size()));

    pairing_choice choice = eval_pairing_choices(s1, s2, score, p1, p2);
    priQueue.push(choice.non_gap);
    priQueue.push(choice.gap_s1);
    priQueue.push(choice.gap_s2);

    // TODO: keep track of cur score, if t_max < cur score then bail

    while (!priQueue.empty())
    {
        if (!cur_pairing.is_set()) {
            cur_pairing = priQueue.top();
            priQueue.pop();

            if (cur_pairing.ft.max <= best_score) {
                cur_pairing = pairing{}; // reset so we can pop
                continue;
            }

            cur_pairings.clear();
            base_offset = cur_pairing.s1offset;
        }

        score += pairing_score(cur_pairing.type);
        if (cur_pairing.s1offset + 1 == s1.size() || cur_pairing.s2_offset + 1 == s2.size())
        {
            best_score = cur_pairing.ft.max;
            for (auto it = cur_pairings.begin(); it != cur_pairings.end(); ++it)
            {
                best_pairings[base_offset] = *it;
                ++base_offset;
            }

            // Add gaps where needed
            if (cur_pairing.s1offset + 1 != s1.size())
            {
                for (auto i = cur_pairing.s1offset; i < s1.size(); ++i, ++base_offset)
                {
                    //best_pairings[best_offset] = final_pairing{i, cur_pairing.s2_offset, cur_pairing.type};
                }
            } else
            {
                for (auto i = cur_pairing.s2_offset; i < s2.size(); ++i, ++base_offset){
                    //best_pairings[best_offset] = final_pairing{cur_pairing.s1offset, i, Gap};
                }
            }

            cur_pairing = pairing{}; // reset so we can pop
            continue;
        }

        // TODO: commit current pairing once we have child pairing that is better then existing
        // TODO: compare current best pairing with existing best pairing, if its worse we need to back
        // track
        choice = eval_pairing_choices(s1, s2, score, p1, p2);
        if (choice.non_gap.ft.max > choice.gap_s1.ft.max && choice.non_gap.ft.max > choice.gap_s2.ft.max)
        {
            pairing_key key{choice.non_gap.s1offset, choice.non_gap.s2_offset};
            auto existing = knownPairings.find(key);
            if (existing != knownPairings.end()) {
                if (choice.non_gap.ft.max <= existing->second) {
                    cur_pairing = pairing{}; // reset so we can pop
                    continue;
                }
            }

            knownPairings[key] = choice.non_gap.ft.max;
            //cur_pairings.emplace_back(choice.non_gap.s1offset, choice.non_gap.s2_offsettype);
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

    int64_t s1offset = 0;
    int64_t s2_offset = 0;
    int64_t s1mut_start = -1;
    int64_t s2_mut_start = -1;
    double total_muts = 0;
    vector<mutation> muts;

    for (int64_t i=0; i < best_pairings.size(); ++i)
    {
        auto p = best_pairings[i];

        if (p.s1 == p.s2)
        {
            if (s1mut_start != -1)
            {
                location h1{s1mut_start, s1offset - s1mut_start};
                location h2{s2_mut_start, s2_offset - s2_mut_start};

                muts.emplace_back(h1, h2);
                s1mut_start = -1;
                total_muts += h1.length;
            }

            s1offset++;
            s2_offset++;
        } else if (s1mut_start != -1)
        {
            s1mut_start = s1offset;
            s2_mut_start = s2_offset;
        } else {
            s1offset++;

            if (p.s2 != static_cast<byte>(0)) // mismatch detection
            {
                s2_offset++;
            }
        }
    }

    return alignment_result{std::move(muts), total_muts / best_pairings.size(), ""};
}

} // dna
