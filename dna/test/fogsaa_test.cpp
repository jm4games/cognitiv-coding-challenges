#include "catch.hpp"
#include "fogsaa.hpp"
#include "fake_stream.hpp"

using namespace dna;

TEST_CASE("fogsaa: Given empty sequences no mutations should exist")
{
    fake_stream s1("", 512);
    fake_stream s2("", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 0);
    REQUIRE(res.simularity_score == 1);
}

TEST_CASE("fogsaa: Given an empty sequence paired with non empty, an error should occur")
{
    fake_stream s1("AA", 512);
    fake_stream s2("", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.error.size() > 0);
}

TEST_CASE("fogsaa: Given two stands with single same base, no mutations should eixst")
{
    fake_stream s1("A", 512);
    fake_stream s2("A", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 0);
    REQUIRE(res.simularity_score == 1);
}
