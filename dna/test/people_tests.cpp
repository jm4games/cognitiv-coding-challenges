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
        T& a_;
        T& b_;

    public:
        align_closure(sequence_aligner<T>& aligner, T& a, T& b)
            : aligner_(aligner), a_(a), b_(b)
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
        auto closure = new align_closure(aligner_, a, b);
        return pool_.enqueue(align_main, closure);
    }
};

SCENARIO("Testing dna between people") {
    GIVEN("Two People")
    {
        WHEN("BLAH") {
            //fake_person bob = std::move(fake_person_factory::new_person_with_dup_chromos());
            //fake_person alice = std::move(fake_person_factory::new_person_with_dup_chromos());

            //thread_pool pool(std::thread::hardware_concurrency());

            //fogsaa_aligner<fake_stream> fogsaa;
            //threaded_alignment_forker<fake_stream> forker(pool, fogsaa);
            //pairwise_aligner<fake_stream> aligner(pool, forker);
            //auto res = aligner.analyze_people_async(bob, alice);
        }
    }
}
