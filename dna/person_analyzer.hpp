#pragma once

#include "person.hpp"
#include "sequence_aligner.hpp"
#include "thread_pool.hpp"

#include <future>

namespace dna {

template<typename T>
concept bool ChromesomeSimliar = requires(T a, T b) {
    a.chromosomes() == b.chromosomes();
};

template<typename T>
concept bool IsPerson = requires(T a) {
    a.chromosomes() == 23;
};

template<HelixStream T>
class AlignmentForker {
public:
    virtual std::future<alignment_result> spawn_alignment(T& a, T& b) = 0;
};

template <HelixStream T>
class pairwise_aligner {
    AlignmentForker<T>& forker_;
    thread_pool& pool_;

public:
    explicit pairwise_aligner(thread_pool& pool, AlignmentForker<T>& forker)
        : forker_(forker), pool_(pool) {};

    // TODO: Use future build in error?
    // TODO: use constraints to enforce person chromo count and that both people have same count.
    template <typename P> requires (Person<P> && IsPerson<P>)
    std::future<std::vector<alignment_result>> analyze_people_async(P& p1, P& p2){
        std::vector<std::future<alignment_result>> futures;
        futures.reserve(p1.chromosomes());

        for (std::size_t i = 0; i < p1.chromosomes(); ++i) {
            //auto my_future = forker_.spawn_alignment(p1.chromosome(i), p2.chromosome(i));
            //futures.push_back(std::move(my_future));
        }

        auto wait_for_all = [my_futures = std::move(futures)] {
            std::vector<std::future<alignment_result>> futures2 = std::move(my_futures);
            std::vector<alignment_result> results;
            results.reserve(futures2.size());

            for (size_t i = 0; i < futures2.size(); ++i) {
                results[i] = futures2[i].get();
            }

            return results;
        };

        return pool_.enqueue(wait_for_all);
    }
};

} // dna
