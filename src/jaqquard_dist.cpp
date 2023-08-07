#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <raptor/argument_parsing/search_arguments.hpp>
#include <raptor/adjust_seed.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/search/do_parallel.hpp>
#include <raptor/search/load_index.hpp>
#include <raptor/search/sync_out.hpp>
#include <raptor/adjust_seed.hpp>

#include <robin_hood.h>

#include "jaqquard_dist.hpp"
#include "options.hpp"

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4; // instead of dna5
};

std::vector<uint64_t> compute_sizes(std::vector<std::string> const & filenames,
                                    smash_options const & options)
{
    std::vector<uint64_t> sizes{};
    robin_hood::unordered_set<uint64_t> hashes{};

    auto get_kmers = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{options.kmer_size}},
                                                   seqan3::window_size{options.kmer_size},
                                                   seqan3::seed{raptor::adjust_seed(options.kmer_size)});

    for (auto const & filename : filenames)
    {
        hashes.clear();
        for (auto && rec : seqan3::sequence_file_input<my_traits>{filename})
            for (auto const hash : rec.sequence() | get_kmers)
                hashes.insert(hash);

        sizes.push_back(hashes.size());
        return sizes;
    }

    return sizes;
}

std::vector<std::string> get_index_filenames(raptor::raptor_index<raptor::index_structure::hibf> const & index)
{
    std::vector<std::string> filenames;
    for (auto const & user_bin_paths : index.bin_path())
    {
        if (user_bin_paths.size() != 1)
            throw std::runtime_error{"No multi file user bins allowed yet;"};
        filenames.push_back(user_bin_paths.front());
    }
    return filenames;
}

void jaqquard_dist(smash_options const & options)
{
    auto index = raptor::raptor_index<raptor::index_structure::hibf>{};

    auto get_kmers = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{options.kmer_size}},
                                                   seqan3::window_size{options.kmer_size},
                                                   seqan3::seed{raptor::adjust_seed(options.kmer_size)});

    raptor::search_arguments arguments{.index_file = options.index_file,
                                       .out_file = options.output_file};

    load_index(index, arguments);

    std::vector<std::string> const index_filenames = get_index_filenames(index);

    std::cerr << "Computing index user bin sizes..." << std::endl;
    std::vector<uint64_t> const index_filename_sizes = compute_sizes(index_filenames, options);

    // raptor::sync_out synced_out_mash{options.output_file.string() + ".mash"};
    raptor::sync_out synced_out{arguments};

    std::cerr << "Writing header line..." << std::endl;
    { // write header line
        std::string line{"#filenames"};
        for (auto const & filename : index_filenames)
        {
            line += '\t';
            line += filename;
        }
        line += '\n';
        synced_out.write(line);
    }

    std::cerr << "Computing distances..." << std::endl;

    // auto worker = [&](size_t const start, size_t const end)
    // {
        auto counter = index.template counting_agent<uint32_t>();

        std::string result_string{};

        for (auto const & filename : options.files /* | seqan3::views::slice(start, end) */)
        {
            result_string.clear();
            result_string += filename;

            robin_hood::unordered_set<uint64_t> hashes{};

            for (auto && rec : seqan3::sequence_file_input<my_traits>{filename})
                for (auto const hash : rec.sequence() | get_kmers)
                    hashes.insert(hash);

            // For all hashes computed for current `filename` count their occurence for each user bin in the HIBF
            auto & result = counter.bulk_count(hashes);

            for (size_t i = 0; i < result.size(); ++i)
            {
                /* A intersect B    =                  result[i]               // #shared-hashes
                 * -------------    =   ------------------------------------
                 *   A union B      =   index_filename_sizes[i] + hashes.size() - result[i]  // #hashes-A + #hashes-B - #shared-hashes
                 */

                auto const dist = static_cast<double>(result[i]) / (index_filename_sizes[i] + hashes.size() - result[i]) - options.fpr;

                result_string += '\t';
                result_string += std::to_string(dist);

                // synced_out_mash.write(filename + "\t" + index.bin_path()[i][0] + "\t" +
                //                       std::to_string(result[i]) + "/" + std::to_string(hashes.size()) + "\t" +
                //                       std::to_string(result[i]) + "/" + std::to_string(index_filename_sizes[i]) + "\n");
            }

            result_string += '\n';
            synced_out.write(result_string);
        }
    // };

    // raptor::do_parallel(worker, options.files.size(), options.threads);
}
