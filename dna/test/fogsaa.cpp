#include "fogsaa.hpp"
#include <algorithm>
#include <limits>
#include <cstdint>
#include <memory>
#include <queue>

using namespace std;

namespace dna
{

const std::byte telemere_seq[]{byte{0x3}, byte{0x3}, byte{0x0}, byte{0x2}, byte{0x2}, byte{0x2}};

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
    GapS1,
    GapS2,
};

static constexpr int64_t pairing_score(pairing_type type)
{
    switch(type)
    {
        case Match:
            return MatchScore;
        case MisMatch:
            return MisMatchScore;
        case GapS1:
        case GapS2:
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
    int64_t pairing_offset = -1;
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

class byte_aligner
{
    using score_cache = unordered_map<pairing_key, int64_t, hash_pairing_key>;
    const vector<byte>& s1_;
    const vector<byte>& s2_;

    unique_ptr<final_pairing[]> best_pairings_;
    unique_ptr<final_pairing[]> cur_pairings_;
    int64_t best_pairings_len_ = 0;

    pairing_choice eval_pairing_choices(const int64_t score, int64_t cur_offset, int64_t p1, int64_t p2) const
    {
        int64_t s1_sz = s1_.size() - BASE_S_OFFSET, s2_sz = s2_.size() - BASE_S_OFFSET;
        pairing_choice result;

        // compare
        int64_t p1n = p1 + 1, p2n = p2 + 1;
        result.non_gap.type = s1_[p1n] == s2_[p2n] ? Match : MisMatch;

        int64_t x1 = s2_sz - p2n, x2 = s1_sz - p1n;
        int64_t score_nxt = score + (result.non_gap.type == Match ? MatchScore : MisMatchScore);
        result.non_gap.score = score_nxt;
        result.non_gap.s1_offset = p1n;
        result.non_gap.s2_offset = p2n;
        result.non_gap.ft = fitness_score{score_nxt + fs_min(x1, x2),  score_nxt + fs_max(x1, x2)};
        result.non_gap.pairing_offset = cur_offset + 1;

        score_nxt = score + GapPenalty;

        // gap s1
        x1 = s2_sz - p2n, x2 = s1_sz - p1;
        result.gap_s1.score = score_nxt;
        result.gap_s1.s1_offset = p1;
        result.gap_s1.s2_offset = p2n;
        result.gap_s1.ft = fitness_score{score_nxt + fs_min(x1, x2), score_nxt + fs_max(x1, x2)};
        result.gap_s1.type = GapS1;
        result.gap_s1.pairing_offset = cur_offset + 1;

        // gap s2
        x1 = s2_sz - p2, x2 = s1_sz - p1n;
         result.gap_s2.score = score_nxt;
        result.gap_s2.s1_offset = p1n;
        result.gap_s2.s2_offset = p2;
        result.gap_s2.ft = fitness_score{score_nxt + fs_min(x1, x2), score_nxt + fs_max(x1, x2)};
        result.gap_s2.type = GapS2;
        result.gap_s2.pairing_offset = cur_offset + 1;

        return result;
    }

    alignment_result get_alignment()
    {
        const int64_t s_offset = static_cast<int64_t>(BASE_S_OFFSET);
        int64_t s1_offset = s_offset;
        int64_t s2_offset = s_offset;
        int64_t s1_mut_start = -1;
        int64_t s2_mut_start = -1;
        double total_muts = 0;
        vector<mutation> muts;

        for (int64_t i= BASE_S_OFFSET; i < best_pairings_len_; ++i)
        {
            auto p = best_pairings_[i];

            if (p.s1 == p.s2)
            {
                if (s1_mut_start != -1) // mutation ended, save it
                {
                    location h1{s1_mut_start - s_offset, s1_offset - s1_mut_start};
                    location h2{s2_mut_start - s_offset, s2_offset - s2_mut_start};

                    muts.emplace_back(h1, h2);
                    s1_mut_start = -1;
                    total_muts += h1.length;
                }
            } else if (s1_mut_start == -1) // try to start mutation
            {
                s1_mut_start = s1_offset;
                s2_mut_start = s2_offset;
            }

            if (p.s1 != gapByte)
                s1_offset++;

            if (p.s2 != gapByte)
                s2_offset++;
        }

        if (s1_mut_start != -1) // mutation ended, save it
        {
            location h1{s1_mut_start - s_offset, s1_offset - s1_mut_start};
            location h2{s2_mut_start - s_offset, s2_offset - s2_mut_start};

            muts.emplace_back(h1, h2);
            total_muts += h1.length;
        }

        // TODO: Better way to score?
        return alignment_result{
            std::move(muts),
            1 - (total_muts / (max(s1_.size(), s2_.size()) - BASE_S_OFFSET)),
            ""
        };
    }

    bool is_candidate(const score_cache& cache, const pairing& pairing, int64_t best_score)
    {
        if (pairing.ft.max < best_score)
            return false;

        pairing_key key = pairing_key{pairing.s1_offset, pairing.s2_offset};
        auto existing = cache.find(key);
        return existing == cache.end() || existing->second < pairing.ft.max;
    }

public:
    byte_aligner(const vector<byte>& s1, const vector<byte>& s2)
        : s1_(s1), s2_(s2), best_pairings_()
    {
        size_t size = max(s1.size(), s2.size()) + (max(s1.size(), s2.size())/ 2); // div 2 accounts for overflow
        best_pairings_ = move(make_unique<final_pairing[]>(size));
        cur_pairings_ = move(make_unique<final_pairing[]>(size));
    }

    alignment_result run_alignment() {
        int64_t best_score = numeric_limits<int64_t>::min();
        int64_t best_min = numeric_limits<int64_t>::min();
        unordered_map<pairing_key, int64_t, hash_pairing_key> best_fit_scores;
        // TODO: replace queue vector with ordered map?
        priority_queue<pairing, vector<pairing>, pairing> pri_queue;
        pairing cur_pairing;

        pairing_choice choice = eval_pairing_choices(0, 0, 0, 0);
        pri_queue.push(choice.non_gap);
        pri_queue.push(choice.gap_s1);
        pri_queue.push(choice.gap_s2);

        while (!pri_queue.empty())
        {
            cur_pairing = pri_queue.top();
            pri_queue.pop();

            if (cur_pairing.ft.max <= best_score)
                break; // we are done, top of queue max can't beat best score

            int64_t base_offset = cur_pairing.pairing_offset;
            bool has_candidate = true;
            int64_t cur_len = base_offset;
            while (has_candidate)
            {
                ++cur_len;

                pairing_key key {cur_pairing.s1_offset, cur_pairing.s2_offset};
                best_fit_scores[key] = cur_pairing.ft.max;
                switch (cur_pairing.type)
                {
                    case Match:
                    case MisMatch:
                      cur_pairings_[cur_pairing.pairing_offset] =
                          final_pairing{s1_[cur_pairing.s1_offset], s2_[cur_pairing.s2_offset]};
                    break;
                    case GapS2:
                      cur_pairings_[cur_pairing.pairing_offset] =
                          final_pairing{s1_[cur_pairing.s1_offset], gapByte};
                    break;
                    default:
                      cur_pairings_[cur_pairing.pairing_offset] =
                          final_pairing{gapByte, s2_[cur_pairing.s2_offset]};
                    break;
                }

                if (cur_pairing.s1_offset + 1 == s1_.size() || cur_pairing.s2_offset +1 == s2_.size())
                {
                    best_score = cur_pairing.score;
                    best_min = max(cur_pairing.ft.min, best_min);

                    for (int64_t i = base_offset; i < cur_len; ++i)
                        best_pairings_[i] = cur_pairings_[i];

                    cur_pairing.s1_offset++;
                    for (int64_t i=cur_pairing.s1_offset; i< s1_.size(); ++i)
                    {
                        best_pairings_[cur_len] = final_pairing{s1_[i], gapByte};
                        ++cur_len;
                        best_score -= GapPenalty;
                    }

                    cur_pairing.s2_offset++;
                    for (int64_t i=cur_pairing.s2_offset; i< s2_.size(); ++i)
                    {
                        best_pairings_[cur_len] = final_pairing{gapByte, s2_[i]};
                        ++cur_len;
                        best_score -= GapPenalty;
                    }

                    best_pairings_len_ = cur_len;
                    break;
                }

                choice = eval_pairing_choices(
                        cur_pairing.score,
                        cur_pairing.pairing_offset,
                        cur_pairing.s1_offset,
                        cur_pairing.s2_offset);
                has_candidate = is_candidate(best_fit_scores, choice.non_gap, best_score);
                if (has_candidate)
                    cur_pairing = choice.non_gap;

                // TODO: make helper function
                if (is_candidate(best_fit_scores, choice.gap_s2, best_score))
                {
                    if (!has_candidate)
                    {
                        cur_pairing = choice.gap_s2;
                    }
                    else if (cur_pairing.ft.max < choice.gap_s2.ft.max)
                    {
                        pri_queue.push(cur_pairing);
                        cur_pairing = choice.gap_s2;
                    } else if (!(  choice.gap_s2.ft.max < cur_pairing.ft.min
                                || choice.gap_s2.ft.max < best_min))
                    {
                        pri_queue.push(choice.gap_s2);
                    }
                    has_candidate = true;
                }

                if (is_candidate(best_fit_scores, choice.gap_s1, best_score))
                {
                    if (!has_candidate)
                    {
                        cur_pairing = choice.gap_s1;
                    }
                    else if (cur_pairing.ft.max < choice.gap_s1.ft.max)
                    {
                        pri_queue.push(cur_pairing);
                        cur_pairing = choice.gap_s1;
                    } else if (!(  choice.gap_s1.ft.max < cur_pairing.ft.min
                                || choice.gap_s1.ft.max < best_min))
                    {
                        pri_queue.push(choice.gap_s1);
                    }
                    has_candidate = true;
                }
            }
        }

        return get_alignment();
    }
};

alignment_result fogsaa::align_bytes(const vector<byte>& s1, const vector<byte>& s2)
{
    byte_aligner aligner(s1, s2);
    return aligner.run_alignment();
}

} // dna
