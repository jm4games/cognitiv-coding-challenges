#include "catch.hpp"
#include "fake_person_factory.hpp"
#include <person.hpp>
#include <iostream>

template<dna::Person P>
class person_tester
{
	std::reference_wrapper<std::ostream> writer_;
	std::reference_wrapper<P> person_;
public:
	person_tester(std::ostream& writer, P& person) :
			writer_(writer),
			person_(person)
	{ }

	void dump_chromosomes()
	{
		for (std::size_t i = 0; i < person_.get().chromosomes(); ++i)
		{
			writer_.get() << "Chromosome " << (i + 1) << ": ";
			auto chrom = person_.get().chromosome(i);

			while (true)
			{
				auto buf = chrom.read();
				if (buf.size() == 0)
					break;
				writer_.get() << buf;
			}

			writer_.get() << "\n";
		}

		writer_.get().flush();
	}
};

TEST_CASE("Fake stream is suitable for testing", "[stream]")
{
	auto data = data::fake();
	fake_stream stream(std::move(data), 128);

	REQUIRE(stream.size() == 1020);
	for (int i = 0; i < 7; i++)
	{
		auto seq = stream.read();
		REQUIRE(seq.size() == 512);
		INFO(seq);
	}

	auto endseq = stream.read();
	REQUIRE(endseq.size() == 496);
	INFO(endseq);

	stream.seek(1018);
	endseq = stream.read();
	REQUIRE(endseq.size() == 8);
	REQUIRE(endseq[0] == dna::G);
	REQUIRE(endseq[1] == dna::G);
	REQUIRE(endseq[2] == dna::G);
	REQUIRE(endseq[3] == dna::T);
	REQUIRE(endseq[4] == dna::G);
	REQUIRE(endseq[5] == dna::T);
	REQUIRE(endseq[6] == dna::C);
	REQUIRE(endseq[7] == dna::C);
}

//TEST_CASE("Fake person fulfills Person concept", "[stream]")
//{
	//fake_person person = std::move(fake_person_factory::new_person_with_dup_chromos());
	//person_tester<fake_person> tester(std::cout, person);
	//tester.dump_chromosomes();
//}

