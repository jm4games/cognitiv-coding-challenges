#include "catch.hpp"
#include "fake_person_factory.hpp"
#include <person.hpp>
#include <person_analyzer.hpp>

using namespace dna;

SCENARIO("Testing dna between people") {
    GIVEN("Two People")
    {
        WHEN("BLAH") {
            fake_person bob = std::move(fake_person_factory::new_person_with_dup_chromos());
            fake_person alice = std::move(fake_person_factory::new_person_with_dup_chromos());

            analyzer my_analyzer;
            //auto res = my_analyzer.analyze_people_async(bob, alice);
        }
    }
}
