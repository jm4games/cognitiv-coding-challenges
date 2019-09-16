#include "catch.hpp"
#include "fogsaa.hpp"
#include "fake_stream.hpp"

using namespace dna;

TEST_CASE("Given empty sequences no mutations should exist")
{
    fake_stream s1("", 512);
    fake_stream s2("", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 0);
}
