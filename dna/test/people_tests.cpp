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

public:
    explicit threaded_alignment_forker(thread_pool& pool, sequence_aligner<T>& aligner)
        : pool_(pool), aligner_(aligner) {}

    future<alignment_result> spawn_alignment(T& a, T& b) override {
        return pool_.enqueue(
                [](sequence_aligner<T>& aligner, T& c, T& d) { return aligner.align(c, d); },
                aligner_,
                a,
                b);
    }
};

SCENARIO("Testing dna between people") {
    GIVEN("Two People")
    {
        WHEN("BLAH") {
            fake_person bob = std::move(fake_person_factory::new_person_with_dup_chromos());
            fake_person alice = std::move(fake_person_factory::new_person_with_dup_chromos());

            thread_pool pool(std::thread::hardware_concurrency());

            fogsaa<fake_stream> fogsaa;
            threaded_alignment_forker<fake_stream> forker(pool, fogsaa);
            pairwise_aligner<fake_stream> aligner(pool, forker);
            auto res = aligner.analyze_people_async(bob, alice);
        }
    }
}
