#pragma once

#include "person.hpp"
#include "sequence_aligner.hpp"
#include "thread_pool.hpp"
#include <future>

namespace dna {

const std::string ChromoMismatchMFErr = "Cannot match an X with a Y chromosome";

template<typename T>
concept bool ChromesomeSimliar = requires(T a, T b) {
    a.chromosomes() == b.chromosomes();
};

template<typename T>
concept bool IsPerson = requires(T a) {
    // TODO: Fix ME, probably need to use value type.
    a.chromosomes() == 23;
};

template<HelixStream T>
class AlignmentForker {
public:
    virtual std::future<alignment_result> spawn_alignment(T&& a, T&& b) = 0;
};

// Analyize two pairs of dna
template <HelixStream T>
class pairwise_aligner {
    using future_results = std::vector<std::future<alignment_result>>;

    AlignmentForker<T>& forker_;
    thread_pool& pool_;

    class futures_closure
    {
        future_results futures_;
    public:
        explicit futures_closure(future_results&& futures) : futures_(std::move(futures)) {}

        std::vector<alignment_result> wait_for_all()
        {
            std::vector<alignment_result> results;
            results.reserve(futures_.size());

            for (size_t i = 0; i < futures_.size(); ++i) {
                results.push_back(std::move(futures_[i].get()));
            }

            return results;
        }
    };

public:
    explicit pairwise_aligner(thread_pool& pool, AlignmentForker<T>& forker)
        : forker_(forker), pool_(pool) {};

    // TODO: Use future build in error?
    // TODO: use constraints to enforce person chromo count of 23
    template <typename P> requires (Person<P> && IsPerson<P>)
    std::future<std::vector<alignment_result>> analyze_people_async(P& p1, P& p2){
        future_results futures;
        futures.reserve(p1.chromosomes());

        for (std::size_t i = 0; i < p1.chromosomes(); ++i)
        {
            T h1 = std::move(p1.chromosome(i));
            T h2 = std::move(p2.chromosome(i));

            if (i == 22)
            {
                // Detect X/Y chromosome mismatch. A Y chromosome has ~57 million bp and an X
                // chromosome has ~156 million bp. If one chromesome is less then 60% the size
                // of other we can safely say we have a identified an X/Y (or a really corrupt strand).
                if (static_cast<double>(std::min(h1.size(), h2.size())) / std::max(h1.size(), h2.size()) < .6)
                {
                    std::promise<alignment_result> p;
                    alignment_result result;
                    result.error = std::string(ChromoMismatchMFErr);
                    p.set_value(result);
                    futures.push_back(std::move(p.get_future()));
                    break;
                }
            }

            std::future<alignment_result> my_future = forker_.spawn_alignment(
                    std::move(h1), std::move(h2));
            futures.push_back(std::move(my_future));
        }

        auto wait_main = [](futures_closure* closure)
                {
                    auto results = std::move(closure->wait_for_all());
                    delete closure;
                    return results;
                };
        return pool_.enqueue(wait_main, new futures_closure(std::move(futures)));
    }
};

} // dna
