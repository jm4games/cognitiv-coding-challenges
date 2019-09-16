#pragma once

#include "sequence_buffer.hpp"

namespace dna
{

template<typename T>
concept bool HelixStream = requires(T a) {
	{ a.seek(1000L) };
	{ a.read() } -> sequence_buffer<ByteBuffer>
	{ a.size() } -> std::size_t;
};

template<typename T>
concept bool Person = requires(T a) {
	{ a.chromosome(1) } -> HelixStream
	{ a.chromosomes() } -> std::size_t
};

} // dna

