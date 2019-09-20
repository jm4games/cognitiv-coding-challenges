#include "catch.hpp"
#include "fogsaa.hpp"
#include "fake_stream.hpp"

using namespace dna;

//TEST_CASE("fogsaa: Given empty sequences no mutations should exist")
//{
    //fake_stream s1("", 512);
    //fake_stream s2("", 512);
    //fogsaa aligner;

    //alignment_result res = aligner.align(s1, s2);
    //REQUIRE(res.mutations.size() == 0);
    //REQUIRE(res.simularity_score == 1);
//}

//TEST_CASE("fogsaa: Given an empty sequence paired with non empty, an error should occur")
//{
    //fake_stream s1("AA", 512);
    //fake_stream s2("", 512);
    //fogsaa aligner;

    //alignment_result res = aligner.align(s1, s2);
    //REQUIRE(res.error.size() > 0);
//}

//TEST_CASE("fogsaa: Given two stands with same bases, no mutations should exist")
//{
    //fake_stream s1("TTTT", 512);
    //fake_stream s2("TTTT", 512);
    //fogsaa aligner;

    //alignment_result res = aligner.align(s1, s2);
    //REQUIRE(res.mutations.size() == 0);
    //REQUIRE(res.simularity_score == 1);
//}

//TEST_CASE("fogsaa: Given two stands with single mutation at end, mutations should exist")
//{
    //fake_stream s1("CCCC", 512);
    //fake_stream s2("CCCT", 512);
    //fogsaa aligner;

    //alignment_result res = aligner.align(s1, s2);
    //REQUIRE(res.mutations.size() == 1);
    //REQUIRE(res.simularity_score == 0.75);
//}

//TEST_CASE("fogsaa: Given two stands with single mutation at begining, mutations should exist")
//{
    //fake_stream s1("TCCC", 512);
    //fake_stream s2("CCCC", 512);
    //fogsaa aligner;

    //alignment_result res = aligner.align(s1, s2);
    //REQUIRE(res.mutations.size() == 1);
    //REQUIRE(res.simularity_score == 0.75);
//}

//TEST_CASE("fogsaa: Given two not similar stands, mutations should exist")
//{
    //fake_stream s1("CCCC", 512);
    //fake_stream s2("AAAA", 512);
    //fogsaa aligner;

    //alignment_result res = aligner.align(s1, s2);
    //REQUIRE(res.mutations.size() == 1);
    //REQUIRE(res.simularity_score == 0);
//}

TEST_CASE("fogsaa: Given two ....")
{
    fake_stream s1("ACGGTTGCCCTT", 512);
    fake_stream s2("AGGGTCCC", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 3);
    for (int i = 0; i < 3; ++i) {
        printf("H1 | Offset: %u, len: %u\n", res.mutations[i].helix1.offset, res.mutations[i].helix1.length);
        printf("H2 | Offset: %u, len: %u\n\n", res.mutations[i].helix2.offset, res.mutations[i].helix2.length);
    }
    REQUIRE(res.simularity_score == 5.0/8.0);
}
