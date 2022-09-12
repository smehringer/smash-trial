#pragma once

#include <vector>
#include <queue>

#include <seqan3/search/views/kmer_hash.hpp>

template <typename T>
struct my_priority_queue : public std::priority_queue<T>
{
    auto get_underlying_container()
    {
        return this->c;
    }
};


template <typename range_t>
void init_sketch(range_t && input, uint8_t const kmer_size, uint32_t const sketch_size, my_priority_queue<uint64_t> & sketch)
{
    auto hashes = input | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmer_size}});
    auto hashes_it = hashes.begin();

    // initialise sketch with the first hashes
    for (size_t i = 0; i < sketch_size && i < hashes.size(); ++i, ++hashes_it)
        sketch.push(*hashes_it);

    while (hashes_it != hashes.end())
    {
        if (*hashes_it < sketch.top())
        {
            sketch.pop();
            sketch.push(*hashes_it);
        }
    }
}

template <typename range_t>
void add_to_sketch(range_t && input, uint8_t const kmer_size, my_priority_queue<uint64_t> & sketch)
{
    auto hashes = input | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmer_size}});
    auto hashes_it = hashes.begin();

    while (hashes_it != hashes.end())
    {
        if (*hashes_it < sketch.top())
        {
            sketch.pop();
            sketch.push(*hashes_it);
        }
    }
}

// can be used if only hashing a single sequence
template <typename range_t>
std::vector<uint64_t> sketch_min_hash(range_t && input, uint8_t const kmer_size, uint32_t const sketch_size)
{
    my_priority_queue<uint64_t> sketch;

    auto hashes = input | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmer_size}});
    auto hashes_it = hashes.begin();

    // initialise sketch with the first hashes
    for (size_t i = 0; i < sketch_size && i < hashes.size(); ++i, ++hashes_it)
        sketch.push(*hashes_it);

    while (hashes_it != hashes.end())
    {
        if (*hashes_it < sketch.top())
        {
            sketch.pop();
            sketch.push(*hashes_it);
        }
    }

    return sketch.get_underlying_container();
}
