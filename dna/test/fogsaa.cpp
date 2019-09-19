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
static const byte gapByte = static_cast<byte>(0);

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
    byte s1 = gapByte;
    byte s2 = gapByte;

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
    pairing non_gap;
    //pairing gap_s1;
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

class byte_aligner
{
    using pairing_cache = unordered_map<pairing_key, int64_t, hash_pairing_key>;
    const vector<byte>& s1_;
    const vector<byte>& s2_;

    vector<final_pairing> best_pairings_;
    vector<final_pairing> cur_pairings_;


    // TODO: replace queue vector with ordered map?
    priority_queue<pairing, vector<pairing>, pairing> priQueue_;

    pairing_choice eval_pairing_choices(int64_t score, int64_t p1, int64_t p2)
    {
        int64_t s1sz = s1_.size(), s2_sz = s2_.size();
        pairing_choice result;

        // compare
        int64_t p1n = p1 + 1, p2n = p2 + 1;
        result.non_gap.type = s1_[p1n] == s2_[p2n] ? Match : MisMatch;

        int64_t x1 = s2_sz - p2n, x2 = s1sz - p1n;
        int64_t score_nxt = score + (result.non_gap.type == Match ? MatchScore : MisMatchScore);
        result.non_gap.s1_offset = p1n;
        result.non_gap.s2_offset = p2n;
        result.non_gap.ft = fitness_score{score_nxt + fs_min(x1, x2),  score_nxt + fs_max(x1, x2)};

        score_nxt = score + GapPenalty;

        // gap s1
        //x1 = s2_sz - p2n, x2 = s1sz - min(p1, (int64_t)0);
        //result.gap_s1.s1_offset = p1;
        //result.gap_s1.s2_offset = p2n;
        //result.gap_s1.ft = fitness_score{score_nxt + fs_min(x1, x2), score_nxt + fs_max(x1, x2)};
        //result.gap_s1.type = GapS1;

        // gap s2
        x1 = s2_sz - min(p2, (int64_t)0), x2 = s1sz - p1n;
        result.gap_s2.s1_offset = p1n;
        result.gap_s2.s2_offset = p2;
        result.gap_s2.ft = fitness_score{score_nxt + fs_min(x1, x2), score_nxt + fs_max(x1, x2)};
        result.gap_s2.type = Gap;

        return result;
    }

    void save_current_path(const pairing& cur_pairing, int64_t base_offset)
    {
        for (auto it = cur_pairings_.begin(); it != cur_pairings_.end(); ++it)
        {
            best_pairings_[base_offset] = *it;
            ++base_offset;
        }

        // Add gaps where needed
        if (cur_pairing.s1_offset + 1 != s1_.size())
        {
            for (auto i = cur_pairing.s1_offset; i < s1_.size(); ++i, ++base_offset)
            {
                best_pairings_[base_offset] = final_pairing{s1_[cur_pairing.s1_offset], gapByte};
            }
        } else
        {
            for (auto i = cur_pairing.s2_offset; i < s2_.size(); ++i, ++base_offset){
                best_pairings_[base_offset] = final_pairing{gapByte, s2_[cur_pairing.s2_offset]};
            }
        }
    }

    bool try_insert_pairing(pairing_cache& cache, const pairing& cur_pairing)
    {
        pairing_key key{cur_pairing.s1_offset, cur_pairing.s2_offset};
        auto existing = cache.find(key);
        if (existing != cache.end()) {
            if (cur_pairing.ft.max <= existing->second) {
                return false;
            }
        }

        cache[key] = cur_pairing.ft.max;
        return true;
    }

    alignment_result get_alignment()
    {
        int64_t s1_offset = 0;
        int64_t s2_offset = 0;
        int64_t s1mut_start = -1;
        int64_t s2_mut_start = -1;
        double total_muts = 0;
        vector<mutation> muts;

        for (int64_t i=0; i < best_pairings_.size(); ++i)
        {
            auto p = best_pairings_[i];

            if (p.s1 == p.s2)
            {
                if (s1mut_start != -1)
                {
                    location h1{s1mut_start, s1_offset - s1mut_start};
                    location h2{s2_mut_start, s2_offset - s2_mut_start};

                    muts.emplace_back(h1, h2);
                    s1mut_start = -1;
                    total_muts += h1.length;
                }

                s1_offset++;
                s2_offset++;
            } else if (s1mut_start != -1)
            {
                s1mut_start = s1_offset;
                s2_mut_start = s2_offset;
            } else {
                s1_offset++;

                if (p.s2 != gapByte) // mismatch detection
                {
                    s2_offset++;
                }
            }
        }

        return alignment_result{std::move(muts), total_muts / best_pairings_.size(), ""};
    }

public:
    byte_aligner(const vector<byte>& s1, const vector<byte>& s2)
        : s1_(s1), s2_(s2)
    {
        best_pairings_.reserve(max(s1.size(), s2.size()));
        cur_pairings_.reserve(max(s1.size(), s2.size()));

        pairing_choice choice = eval_pairing_choices(0, -1, -1);
        priQueue_.push(choice.non_gap);
        //priQueue_.push(choice.gap_s1);
        priQueue_.push(choice.gap_s2);
    }

    alignment_result run_alignment() {
        int64_t score = 0;
        int64_t best_score = numeric_limits<int64_t>::min();
        pairing_cache known_pairings;
        unordered_map<pairing_key, int64_t, hash_pairing_key> best_scores;
        pairing cur_pairing;

        while (!priQueue_.empty())
        {
            cur_pairing = priQueue_.top();
            priQueue_.pop();

            if (cur_pairing.ft.max <= best_score) {
                break; // we are done, top of queue max can't beat best score
            }

            // TODO: how do i reset score?
            int64_t base_offset = 0;
            while (true)
            {
                score += pairing_score(cur_pairing.type);

                // end of strand?
                if (cur_pairing.s1_offset + 1 == s1_.size() || cur_pairing.s2_offset + 1 == s2_.size())
                {
                    best_score = score;
                    save_current_path(cur_pairing, base_offset);
                    break;
                }

                pairing_key key{cur_pairing.s1_offset, cur_pairing.s2_offset};
                auto existing = best_scores.find(key);
                if (existing != best_scores.end()) {
                    if (score <= existing->second) {
                        break; // better path exist
                    }
                }

                best_scores[key] = score;
                // TODO: Set cur pairings....

                pairing_choice choice = eval_pairing_choices(
                        score, cur_pairing.s1_offset, cur_pairing.s2_offset);
                pairing other;

                if (choice.non_gap.ft.max >= choice.gap_s2.ft.max) {
                    cur_pairing = choice.non_gap;
                    other = choice.gap_s2;
                } else
                {
                    cur_pairing = choice.gap_s2;
                    other = choice.non_gap;
                }
                // add cur score to pairing
                // we can double add to queue for a position if tmax is better (just check on pop if
                // value has been evaluated alrdy)

                // TODO: I might just needed a visted check on entry to query (independent of tmax)....

                if (  cur_pairing.ft.max <= best_score
                   || !try_insert_pairing(known_pairings, cur_pairing)) {
                    break; // TODO: how do i recover state (score/node) here?
                }

                if (try_insert_pairing(known_pairings, other)) {
                    priQueue_.push(other);
                }
            }
        }

        return get_alignment();
    }
};

alignment_result fogsaa::align_bytes(const vector<byte>& s1, const vector<byte>& s2) const
{
    byte_aligner aligner(s1, s2);
    return aligner.run_alignment();
}

} // dna
