#pragma once

#include "person.hpp"
#include "sequence_aligner.hpp"
#include <string.h>

namespace dna
{

static const size_t BASE_S_OFFSET = 1;

// FOGSAA based on the following whitepaper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3638164/#s1

class fogsaa {

    static const std::byte telemere_seq[6];

    template<HelixStream T>
    static constexpr sequence_buffer_iterator<T>consume_it(sequence_buffer<T>& seq, int64_t count)
    {
        auto it = seq.begin;
        for(;count > 0; --count)
            ++it;

        return it;
    }

    template<HelixStream T>
    static constexpr sequence_buffer_iterator<T> trim_telemere_front(T& helix, sequence_buffer<T>& seq)
    {
        const size_t t_size = sizeof(telemere_seq);
        int64_t chromo_start = 0;
        size_t offset = 0;
        bool find_start = true;

        while (true)
        {
            seq = helix.read();
            if (seq.size() == 0)
                break;

            std::byte* buffer = (&(seq.buffer()[0])) + offset;

            if (find_start)
            {
                // We could support trimming buffers smaller then the t_size, but that makes trimming
                // more complicated. In practice we will never have a buffer smaller then t_size having
                // telemeres so this should not be an issue.
                if (seq.size() < t_size)
                {
                    return seq.begin();
                }

                for (size_t g=0; g < t_size; ++g)
                {
                    if (memcmp(buffer+g, telemere_seq, t_size) == 0)
                    {
                        if (g > 0)
                            find_start = !(memcmp(buffer, telemere_seq + t_size - g, g) == 0);// check for partial match
                        else
                            find_start = false;

                        offset = g + t_size;
                        chromo_start = g + t_size;
                        buffer = buffer + chromo_start;
                        break;
                    }
                }

                if (find_start)
                {
                    return seq.begin(); // no telemeres found at start of strand
                }
            }

            size_t size = seq.size() - offset;
            size_t limit = ((size / t_size) * t_size) + 1;
            for (size_t i = 0; i < limit; i = i + t_size)
            {
                if (memcmp(buffer + i, telemere_seq, std::min(limit - i, t_size)) != 0)
                    return consume_it(seq, chromo_start);
                chromo_start += t_size;
            }

            offset = size - limit - 1; // TODO: do i need -1 here?
        }

        return consume_it(seq, chromo_start);
    }

    template<HelixStream T>
    static void trim_telemere_end(std::vector<std::byte>& vec)
    {

    }

    template<HelixStream T>
    static constexpr void fill_helix_vector(T& helix, std::vector<std::byte>& vec)
    {
        vec.reserve(helix.size() + BASE_S_OFFSET);
        for (size_t i = 0; i < BASE_S_OFFSET; ++i)
            vec.emplace_back(static_cast<std::byte>(0));

        bool trim_start = true;

        while (true)
        {
            auto seq = helix.read();
            if (seq.size() == 0)
                return;

            auto it = seq.begin();
            if (trim_start)
            {
                it = trim_telemere_front(helix, seq);
                trim_start = false;
            }

            auto end = seq.end();
            for (; it != end; ++it)
            {
                vec.emplace_back(static_cast<std::byte>(*it));
            }
        }
    }

    static alignment_result align_bytes(
            const std::vector<std::byte>& s1, const std::vector<std::byte>& s2);
public:

    template<HelixStream T>
    static alignment_result align(T& stream1, T& stream2)
    {
        if (stream1.size() == 0 && stream2.size() == 0)
        {
            alignment_result res;
            res.simularity_score = 1;
            return res;
        }

        if (stream1.size() == 0 || stream2.size() == 0)
        {
            alignment_result res;
            res.error = "A stream did not have data";
            return res;
        }

        std::vector<std::byte> s1;
        std::vector<std::byte> s2;

        // TODO: set longest stream to s1, then flip back as needed
        fill_helix_vector(stream1, s1);
        fill_helix_vector(stream2, s2);

        return align_bytes(s1, s2);
    }
};

template<HelixStream T>
class fogsaa_aligner : public sequence_aligner<T>
{
public:
    fogsaa_aligner() {};

    alignment_result align(T& a, T&b) const override {
        return fogsaa::align(a, b);
    }
};

} // dna
