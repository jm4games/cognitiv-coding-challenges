#include "catch.hpp"
#include "fake_person_factory.hpp"
#include "fogsaa.hpp"
#include <person.hpp>
#include <person_analyzer.hpp>

using namespace dna;
using namespace std;

template<HelixStream T>
class threaded_alignment_forker : public AlignmentForker<T> {
    sequence_aligner<T>& aligner_;
    thread_pool& pool_;

    class align_closure
    {
        sequence_aligner<T>& aligner_;
        T a_;
        T b_;

    public:
        align_closure(sequence_aligner<T>& aligner, T&& a, T&& b)
            : aligner_(aligner), a_(std::move(a)), b_(std::move(b))
        {}

        alignment_result do_align()
        {
            return aligner_.align(a_, b_);
        }
    };

    static alignment_result align_main(align_closure* closure)
    {
        alignment_result res = closure->do_align();
        delete closure;
        return res;
    }

public:
    explicit threaded_alignment_forker(thread_pool& pool, sequence_aligner<T>& aligner)
        : pool_(pool), aligner_(aligner) {}

    future<alignment_result> spawn_alignment(T&& a, T&& b) override {
        auto closure = new align_closure(aligner_, std::move(a), std::move(b));
        return pool_.enqueue(align_main, closure);
    }
};

TEST_CASE("Given two male people with exact dna match, no mutations should be found") {
    fake_person bob = std::move(fake_person_factory::new_person_with_dup_chromos());
    fake_person john = std::move(fake_person_factory::new_person_with_dup_chromos());

    thread_pool pool(std::thread::hardware_concurrency());
    fogsaa_aligner<fake_stream> fogsaa;
    threaded_alignment_forker<fake_stream> forker(pool, fogsaa);
    pairwise_aligner<fake_stream> aligner(pool, forker);
    auto results = aligner.analyze_people_async(bob, john).get();

    REQUIRE(results.size() == 23);

    size_t i = 0;
    for (auto it = results.begin(); it != results.end(); ++it)
    {
        UNSCOPED_INFO("Chromosome " << i);
        alignment_result res = *it;
        REQUIRE(res.mutations.size() == 0);
        REQUIRE(res.similarity_score == 1);
        ++i;
    }
}

TEST_CASE("Given one male and one female, chromosome 23 should have error") {
    fake_person bob = std::move(fake_person_factory::new_person_with_dup_chromos_male());
    fake_person alice = std::move(fake_person_factory::new_person_with_dup_chromos());

    thread_pool pool(std::thread::hardware_concurrency());
    fogsaa_aligner<fake_stream> fogsaa;
    threaded_alignment_forker<fake_stream> forker(pool, fogsaa);
    pairwise_aligner<fake_stream> aligner(pool, forker);
    auto results = aligner.analyze_people_async(bob, alice).get();

    REQUIRE(results.size() == 23);

    size_t i = 0;
    for (auto it = results.begin(); it != results.end(); ++it)
    {
        alignment_result res = *it;
        UNSCOPED_INFO("Chromosome " << i);

        if (i != 22)
        {
            REQUIRE(res.mutations.size() == 0);
            REQUIRE(res.similarity_score == 1);
        } else
        {
            REQUIRE(res.error == dna::ChromoMismatchMFErr);
        }

        ++i;
    }
}
