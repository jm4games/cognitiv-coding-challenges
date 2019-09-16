#include "catch.hpp"
#include "fake_person_factory.hpp"
#include <person.hpp>
#include <person_analyzer.hpp>

using namespace dna;
using namespace std;

template<HelixStream T>
class threaded_alignment_forker : public AlignmentForker<T> {
    SequenceAligner<T>& aligner_;

public:
    explicit threaded_alignment_forker(SequenceAligner<T>& aligner) : aligner_(aligner) {}

    future<alignment_result> spawn_alignment(T& a, T& b) override {
        promise<alignment_result> my_promise;
        // TODO: SPAWN ALIGMENT

        return my_promise.get_future();
    }
};

SCENARIO("Testing dna between people") {
    GIVEN("Two People")
    {
        WHEN("BLAH") {
            fake_person bob = std::move(fake_person_factory::new_person_with_dup_chromos());
            fake_person alice = std::move(fake_person_factory::new_person_with_dup_chromos());

            // TODO: construct fogsaa instance
            threaded_alignment_forker<fake_stream> forker;
            pairwise_aligner<fake_stream> aligner(forker);
            //auto res = my_analyzer.analyze_people_async(bob, alice);
        }
    }
}
