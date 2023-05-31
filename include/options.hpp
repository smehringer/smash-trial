#pragma once

#include <vector>
#include <filesystem>
#include <string>

#include <robin_hood.h>

struct smash_options
{
    std::filesystem::path input_file{};
    std::filesystem::path index_file{};
    std::filesystem::path output_file{};
    uint32_t sketch_size{10000};
    uint8_t kmer_size{32};
    double fpr{0.0};
    uint8_t threads{32};
    bool write_time{true};
    bool no_sketching{false};

    // data
    std::vector<std::string> files;
    robin_hood::unordered_map<std::string, uint64_t> sizes;
};
