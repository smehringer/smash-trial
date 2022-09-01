#pragma once

#include <vector>
#include <filesystem>
#include <string>

struct smash_options
{
    std::filesystem::path input_file{};
    std::filesystem::path index_file{};
    std::filesystem::path output_file{};
    uint32_t sketch_size{10000};
    uint8_t kmer_size{32};
    uint8_t threads{32};
    bool write_time{true};

    // data
    std::vector<std::string> files;
};
