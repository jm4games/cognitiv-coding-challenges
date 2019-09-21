#include "catch.hpp"
#include "fogsaa.hpp"
#include "fake_stream.hpp"

using namespace dna;

void print_mutations(std::vector<mutation>& muts)
{
    for (int i = 0; i < muts.size(); ++i) {
        printf("H1 | Offset: %u, len: %u\n", muts[i].helix1.offset, muts[i].helix1.length);
        printf("H2 | Offset: %u, len: %u\n\n", muts[i].helix2.offset, muts[i].helix2.length);
    }
}

// TODO: Validate mutations

TEST_CASE("Given empty sequences no mutations should exist")
{
    fake_stream s1("", 512);
    fake_stream s2("", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 0);
    REQUIRE(res.simularity_score == 1);
}

TEST_CASE("Given an empty sequence paired with non empty, an error should occur")
{
    fake_stream s1("AA", 512);
    fake_stream s2("", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.error.size() > 0);
}

TEST_CASE("Given two stands with single same base, no mutations should exist")
{
    fake_stream s1("T", 512);
    fake_stream s2("T", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 0);
    REQUIRE(res.simularity_score == 1);
}

TEST_CASE("Given two stands with same bases, no mutations should exist")
{
    fake_stream s1("TTTT", 512);
    fake_stream s2("TTTT", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 0);
    REQUIRE(res.simularity_score == 1);
}

TEST_CASE("Given two stands with single mutation at end, mutations should exist")
{
    fake_stream s1("CCCC", 512);
    fake_stream s2("CCCT", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 1);
    REQUIRE(res.simularity_score == 0.75);
}

TEST_CASE("Given two stands with single mutation at begining, mutations should exist")
{
    fake_stream s1("TCCC", 512);
    fake_stream s2("CCCC", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 1);
    REQUIRE(res.simularity_score == 0.75);
}

TEST_CASE("Given two not similar stands, mutations should exist")
{
    fake_stream s1("CCCC", 512);
    fake_stream s2("AAAA", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 1);
    REQUIRE(res.simularity_score == 0);
}

TEST_CASE("Given two not similar stands (center mutation), 1 mutation should exist")
{
    fake_stream s1("CTTC", 512);
    fake_stream s2("CAAC", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 1);
    REQUIRE(res.simularity_score == 0.50);
}

TEST_CASE("Given two not similar stands (alt mutations), 2 mutations should exist")
{
    fake_stream s1("ATTC", 512);
    fake_stream s2("CTAC", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 2);
    REQUIRE(res.simularity_score == 0.50);
}

TEST_CASE("Given two similar stands where one is offset, 1 mutation should exist")
{
    fake_stream s1("AT", 512);
    fake_stream s2("T", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    REQUIRE(res.mutations.size() == 1);
    REQUIRE(res.simularity_score == 0.50);
}

TEST_CASE("Given two similar stands with one gap, 2 mutation should exist")
{
    fake_stream s1("ACGG", 512);
    fake_stream s2("AGG", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    print_mutations(res.mutations);
    REQUIRE(res.mutations.size() == 1);
    REQUIRE(res.simularity_score == 0.75);
}

TEST_CASE("Given two similar stands with one gap and 1 mismatch (at start), 1 mutation should exist")
{
    fake_stream s1("ACGG", 512);
    fake_stream s2("TGG", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    print_mutations(res.mutations);
    REQUIRE(res.mutations.size() == 1);
    REQUIRE(res.simularity_score == 0.50);
}

TEST_CASE("Given two similar stands with one gap and 1 mismatch (at end), 1 mutation should exist")
{
    fake_stream s1("ACGG", 512);
    fake_stream s2("ACT", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    print_mutations(res.mutations);
    REQUIRE(res.mutations.size() == 1);
    REQUIRE(res.simularity_score == 0.50);
}

TEST_CASE("Given two ....")
{
    fake_stream s1("ACGGTTGC", 512);
    fake_stream s2("AGCGTC", 512);
    fogsaa aligner;

    alignment_result res = aligner.align(s1, s2);
    print_mutations(res.mutations);
    REQUIRE(res.mutations.size() == 3);
    REQUIRE(res.simularity_score == 0.5);
}
