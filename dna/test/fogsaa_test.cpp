#include "catch.hpp"
#include "fogsaa.hpp"
#include "fake_stream.hpp"
#include "fake_person_factory.hpp"

using namespace dna;

void print_mutations(std::vector<mutation>& muts)
{
    for (int i = 0; i < muts.size(); ++i) {
        printf("H1 | Offset: %u, len: %u\n", muts[i].helix1.offset, muts[i].helix1.length);
        printf("H2 | Offset: %u, len: %u\n\n", muts[i].helix2.offset, muts[i].helix2.length);
    }
}

//TEST_CASE("Given empty sequences no mutations should exist")
//{
    //fake_stream s1("", 512);
    //fake_stream s2("", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 0);
    //REQUIRE(res.similarity_score == 1);
//}

//TEST_CASE("Given an empty sequence paired with non empty, an error should occur")
//{
    //fake_stream s1("AA", 512);
    //fake_stream s2("", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.error.size() > 0);
//}

//TEST_CASE("Given two stands with single same base, no mutations should exist")
//{
    //fake_stream s1("T", 512);
    //fake_stream s2("T", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 0);
    //REQUIRE(res.similarity_score == 1);
//}

//TEST_CASE("Given two stands with same bases, no mutations should exist")
//{
    //fake_stream s1("TTTT", 512);
    //fake_stream s2("TTTT", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 0);
    //REQUIRE(res.similarity_score == 1);
//}

//TEST_CASE("Given two stands with single mutation at end, mutations should exist")
//{
    //fake_stream s1("CCCC", 512);
    //fake_stream s2("CCCT", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 1);
    //REQUIRE(res.similarity_score == 0.75);

    //mutation mut{location{3,1}, location{3,1}};
    //REQUIRE(res.mutations[0] == mut);
//}

//TEST_CASE("Given two stands with single mutation at begining, mutations should exist")
//{
    //fake_stream s1("TCCC", 512);
    //fake_stream s2("CCCC", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 1);
    //REQUIRE(res.similarity_score == 0.75);

    //mutation mut{location{0,1}, location{0,1}};
    //REQUIRE(res.mutations[0] == mut);
//}

//TEST_CASE("Given two not similar stands, mutations should exist")
//{
    //fake_stream s1("CCCC", 512);
    //fake_stream s2("AAAA", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 1);
    //REQUIRE(res.similarity_score == 0);

    //mutation mut{location{0,4}, location{0,4}};
    //REQUIRE(res.mutations[0] == mut);
//}

//TEST_CASE("Given two not similar stands (center mutation), 1 mutation should exist")
//{
    //fake_stream s1("CTTC", 512);
    //fake_stream s2("CAAC", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 1);
    //REQUIRE(res.similarity_score == 0.50);

    //mutation mut{location{1,2}, location{1,2}};
    //REQUIRE(res.mutations[0] == mut);
//}

//TEST_CASE("Given two not similar stands (alt mutations), 2 mutations should exist")
//{
    //fake_stream s1("ATTC", 512);
    //fake_stream s2("CTAC", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 2);
    //REQUIRE(res.similarity_score == 0.50);

    //mutation mut{location{0,1}, location{0,1}};
    //REQUIRE(res.mutations[0] == mut);

    //mutation mut2{location{2,1}, location{2,1}};
    //REQUIRE(res.mutations[1] == mut2);
//}

//TEST_CASE("Given two similar stands where one is offset, 1 mutation should exist")
//{
    //fake_stream s1("AT", 512);
    //fake_stream s2("T", 512);

    //alignment_result res = fogsaa::align(s1, s2);

    //REQUIRE(res.mutations.size() == 1);
    //REQUIRE(res.similarity_score == 0.50);

    //mutation mut{location{0,1}, location{0,0}};
    //REQUIRE(res.mutations[0] == mut);
//}

//TEST_CASE("Given two similar stands with one gap, 1 mutation should exist")
//{
    //fake_stream s1("ACGG", 512);
    //fake_stream s2("AGG", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 1);
    //REQUIRE(res.similarity_score == 0.75);

    //mutation mut{location{1,1}, location{1,0}};
    //REQUIRE(res.mutations[0] == mut);
//}

//TEST_CASE("Given two similar stands with one gap and 1 mismatch (at start), 1 mutation should exist")
//{
    //fake_stream s1("ACGG", 512);
    //fake_stream s2("TGG", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 1);
    //REQUIRE(res.similarity_score == 0.50);

    //mutation mut{location{0,2}, location{0,1}};
    //REQUIRE(res.mutations[0] == mut);
//}

//TEST_CASE("Given two similar stands with one gap and 1 mismatch (at end), 1 mutation should exist")
//{
    //fake_stream s1("ACGG", 512);
    //fake_stream s2("ACT", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 1);
    //REQUIRE(res.similarity_score == 0.50);

    //mutation mut{location{2,2}, location{2,1}};
    //REQUIRE(res.mutations[0] == mut);
//}

//TEST_CASE("Given two strands with 2 gaps and 1 mismatch, 3 mutations should exist")
//{
    //fake_stream s1("ACGGTTGC", 512);
    //fake_stream s2("AGCGTC", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 3);
    //REQUIRE(res.similarity_score == 0.5);

    //mutation mut{location{1,1}, location{1,0}};
    //REQUIRE(res.mutations[0] == mut);

    //mutation mut2{location{3,2}, location{2,2}};
    //REQUIRE(res.mutations[1] == mut2);

    //mutation mut3{location{6,1}, location{5,0}};
    //REQUIRE(res.mutations[2] == mut3);
//}

//TEST_CASE("Given two stands with gaps on both strands, 2 mutation should exist")
//{
    //fake_stream s1("AATC", 512);
    //fake_stream s2("ACAT", 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 2);
    //REQUIRE(res.similarity_score == 0.75);

    //mutation mut{location{1,0}, location{1,1}};
    //REQUIRE(res.mutations[0] == mut);

    //mutation mut2{location{3,1}, location{4,0}};
    //REQUIRE(res.mutations[1] == mut2);
//}

//TEST_CASE("Given two stands with same large fake data, no mutations should exist")
//{
    //fake_stream s1(data::fake(), 512);
    //fake_stream s2(data::fake(), 512);

    //alignment_result res = fogsaa::align(s1, s2);
    //REQUIRE(res.mutations.size() == 0);
    //REQUIRE(res.similarity_score == 1.0);
//}
