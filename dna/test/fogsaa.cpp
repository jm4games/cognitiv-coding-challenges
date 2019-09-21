#include "fogsaa.hpp"
#include <algorithm>
#include <limits>
#include <cstdint>
#include <memory>
#include <queue>

using namespace std;

namespace dna
{

static const int64_t MatchScore = 1;
static const int64_t MisMatchScore = -1;
static const int64_t GapPenalty = -2;
static const byte gapByte = static_cast<byte>(8);

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

    final_pairing() {};

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
    int64_t score = numeric_limits<int64_t>::min();
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

    unique_ptr<final_pairing[]> best_pairings_;
    unique_ptr<final_pairing[]> cur_pairings_;

    pairing_choice eval_pairing_choices(const int64_t score, const  int64_t p1, const int64_t p2) const
    {
        int64_t s1_sz = s1_.size(), s2_sz = s2_.size();
        pairing_choice result;

        // compare
        int64_t p1n = p1 + 1, p2n = p2 + 1;
        result.non_gap.type = s1_[p1] == s2_[p2] ? Match : MisMatch;

        int64_t x1 = s2_sz - p2n, x2 = s1_sz - p1n;
        int64_t score_nxt = score + (result.non_gap.type == Match ? MatchScore : MisMatchScore);
        result.non_gap.score = score_nxt;
        result.non_gap.s1_offset = p1;
        result.non_gap.s2_offset = p2;
        result.non_gap.ft = fitness_score{score_nxt + fs_min(x1, x2),  score_nxt + fs_max(x1, x2)};

        score_nxt = score + GapPenalty;

        // gap s1
        //x1 = s2_sz - p2n, x2 = s1_sz - min(p1, (int64_t)0);
        //result.gap_s1.s1_offset = p1;
        //result.gap_s1.s2_offset = p2n;
        //result.gap_s1.ft = fitness_score{score_nxt + fs_min(x1, x2), score_nxt + fs_max(x1, x2)};
        //result.gap_s1.type = GapS1;

        // gap s2
        x1 = s2_sz - p2, x2 = s1_sz - p1n;
         result.gap_s2.score = score_nxt;
        result.gap_s2.s1_offset = p1;
        result.gap_s2.s2_offset = p2;
        result.gap_s2.ft = fitness_score{score_nxt + fs_min(x1, x2), score_nxt + fs_max(x1, x2)};
        result.gap_s2.type = Gap;

        return result;
    }

    bool is_expanded(pairing_cache& cache, const pairing& cur_pairing) const
    {
        pairing_key key{cur_pairing.s1_offset, cur_pairing.s2_offset};
        return cache.find(key) != cache.end();
    }

    alignment_result get_alignment()
    {
        int64_t s1_offset = 0;
        int64_t s2_offset = 0;
        int64_t s1_mut_start = -1;
        int64_t s2_mut_start = -1;
        double total_muts = 0;
        vector<mutation> muts;

        for (int64_t i=0; i < s1_.size(); ++i)
        {
            auto p = best_pairings_[i];

            if (p.s1 == p.s2)
            {
                if (s1_mut_start != -1) // mutation ended, save it
                {
                    location h1{s1_mut_start, s1_offset - s1_mut_start};
                    location h2{s2_mut_start, s2_offset - s2_mut_start};

                    muts.emplace_back(h1, h2);
                    s1_mut_start = -1;
                    total_muts += h1.length;
                }
            } else if (s1_mut_start == -1) // try to start mutation
            {
                s1_mut_start = s1_offset;
                s2_mut_start = s2_offset;
            }

            s1_offset++;

            if (p.s2 != gapByte)
            {
                s2_offset++;
            }
        }

        if (s1_mut_start != -1) // mutation ended, save it
        {
            location h1{s1_mut_start, s1_offset - s1_mut_start};
            location h2{s2_mut_start, s2_offset - s2_mut_start};

            muts.emplace_back(h1, h2);
            total_muts += h1.length;
        }

        return alignment_result{std::move(muts), 1 - (total_muts / s1_.size()), ""};
    }

public:
    byte_aligner(const vector<byte>& s1, const vector<byte>& s2)
        : s1_(s1), s2_(s2), best_pairings_()
    {
        size_t size = max(s1.size(), s2.size());
        best_pairings_ = move(make_unique<final_pairing[]>(size));
        cur_pairings_ = move(make_unique<final_pairing[]>(size));
    }

    alignment_result run_alignment() {
        int64_t best_score = numeric_limits<int64_t>::min();
        unordered_map<pairing_key, int64_t, hash_pairing_key> best_scores;
        // TODO: replace queue vector with ordered map?
        priority_queue<pairing, vector<pairing>, pairing> pri_queue;
        pairing cur_pairing;

        pairing_choice choice = eval_pairing_choices(0, 0, 0);
        pri_queue.push(choice.non_gap);
        //pri_queue.push(choice.gap_s1);
        pri_queue.push(choice.gap_s2);

        while (!pri_queue.empty())
        {
            cur_pairing = pri_queue.top();
            pri_queue.pop();

            if (cur_pairing.ft.max <= best_score)
                break; // we are done, top of queue max can't beat best score

            int64_t base_offset = cur_pairing.s1_offset;
            while (true)
            {
                pairing_key key{cur_pairing.s1_offset, cur_pairing.s2_offset};
                auto existing = best_scores.find(key);
                if (existing != best_scores.end() && cur_pairing.score <= existing->second)
                     break; // better path exist

                best_scores[key] = cur_pairing.score;
                if (cur_pairing.type != Gap)
                {
                  cur_pairings_[cur_pairing.s1_offset] =
                      final_pairing{s1_[cur_pairing.s1_offset], s2_[cur_pairing.s2_offset]};
                } else
                {
                  cur_pairings_[cur_pairing.s1_offset] =
                      final_pairing{s1_[cur_pairing.s1_offset], gapByte};
                }

                // end of strand? remember s2 will always be shortest
                if (cur_pairing.s2_offset + 1 == s2_.size() && cur_pairing.type != Gap)
                {
                    best_score = cur_pairing.score;
                    for (int64_t i = base_offset; i < s1_.size(); ++i)
                        best_pairings_[i] = cur_pairings_[i];
                    break;
                }

                choice = eval_pairing_choices(
                        cur_pairing.score,
                        cur_pairing.s1_offset + 1,
                        cur_pairing.s2_offset + (cur_pairing.type == Gap ? 0 : 1));
                pairing other;

                if (choice.non_gap.ft.max >= choice.gap_s2.ft.max)
                {
                    cur_pairing = choice.non_gap;
                    other = choice.gap_s2;
                } else
                {
                    cur_pairing = choice.gap_s2;
                    other = choice.non_gap;
                }

                key = pairing_key{cur_pairing.s1_offset, cur_pairing.s2_offset};
                existing = best_scores.find(key);

                // TODO: don't enqueue if max is less then best known min
                if (  other.ft.max > best_score && (existing == best_scores.end()
                   || existing->second < other.score))
                {
                    // TODO: do i really want to push other all the time?
                    pri_queue.push(other);
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
