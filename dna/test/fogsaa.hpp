#pragma once

#include "person.hpp"
#include "sequence_aligner.hpp"

namespace dna
{

// FOGSAA based on the following whitepaper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3638164/#s1

template<HelixStream T>
class fogsaa : public sequence_aligner<T> {

    static const std::byte telemere_seq[];

    static constexpr void trim_telemere_front(T& helix, std::vector<std::byte> &vec)
    {
        //int64_t telemere_start = 0; // REM: we have offset of atleast 1 before data, so we can use 0
        //std::byte buffer[sizeof(telemere_seq) * 2];
        //size_t i = 0; 
        //sequence_buffer<T> buf;

        //while (true)
        //{
            //buf = helix.read();
            //if (buf.size() == 0)
                //break;

            //for (;i < sizeof(buffer); ++i)
            //{
                //buffer[i] = static_cast<std::byte>(*it));
            //}
    }

    static constexpr void fill_helix_vector(T& helix, std::vector<std::byte> &vec)
    {
        vec.reserve(helix.size() + BASE_S_OFFSET);
        for (size_t i = 0; i < BASE_S_OFFSET; ++i)
            vec.emplace_back(static_cast<std::byte>(0));

        while (true)
        {
            auto buf = helix.read();
            if (buf.size() == 0)
                return;

            auto end = buf.end();
            for (auto it = buf.begin(); it != end; ++it)
            {
                vec.emplace_back(static_cast<std::byte>(*it));
            }
        }
    }

    alignment_result align_bytes(
            const std::vector<std::byte>& s1, const std::vector<std::byte>& s2) const;
public:
    static const size_t BASE_S_OFFSET;

    fogsaa() {};

    alignment_result align(T& stream1, T& stream2) override
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

} // dna
