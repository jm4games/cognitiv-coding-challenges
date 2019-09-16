#pragma once

#include "sequence_buffer.hpp"
#include <cstdint>
#include <vector>

namespace dna {

// Defines the location of a subsequence in a helix.
struct Location
{
    int64_t offset;
    int64_t size;
};

// Defines the location where two strands are conflicting.
struct Mutation
{
    Location helix1;
    Location helix2;
};

// Describes how similar two helixes are and where they differ (mutations).
struct AlignmentReport
{
    std::vector<Mutation> mutations;
    double simularityScore = 0;
};

//template<template T>
//concept bool SequenceAligner = requires(T a) {
	//{ a.alignWith(T) } -> AlignmentReport
//};

}
