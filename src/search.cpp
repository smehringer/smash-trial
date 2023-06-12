#include <seqan3/search/views/kmer_hash.hpp>

#include <chopper/sketch/hyperloglog.hpp>
#include <chopper/sketch/execute.hpp>
#include <chopper/sketch/estimate_kmer_counts.hpp>

#include <raptor/argument_parsing/search_arguments.hpp>
#include <raptor/adjust_seed.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/search/do_parallel.hpp>
#include <raptor/search/load_index.hpp>
#include <raptor/search/sync_out.hpp>

#include "compute_distance.hpp"
#include "search.hpp"
#include "options.hpp"
#include "sketch.hpp"

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4; // instead of dna5
};

void search(smash_options & options)
{
    auto index = raptor::raptor_index<raptor::index_structure::hibf>{};

    raptor::search_arguments arguments{.index_file = options.index_file,
                                       .out_file = options.output_file};

    raptor::load_index(index, arguments);

    {
        std::vector<std::string> files = options.files;
        std::vector<chopper::sketch::hyperloglog> sketches{};

        chopper::configuration config{.data_file = options.input_file,
                                    .k = options.kmer_size,
                                    .disable_sketch_output = true,
                                    .threads = options.threads};

        for (size_t i = 0; i < index.bin_path().size(); ++i)
        {
            files.push_back(index.bin_path()[i][0]);
            if (index.bin_path()[i].size() > 1)
                throw std::runtime_error{"Multi file user bins not supported yet."};
        }

        chopper::sketch::execute(config, files, sketches);
        std::vector<size_t> kmer_counts{};
        chopper::sketch::estimate_kmer_counts(sketches, kmer_counts);

        for (size_t i = 0; i < files.size(); ++i)
            options.sizes.emplace(files[i], static_cast<uint64_t>(kmer_counts[i]));
    }

    raptor::sync_out synced_out{arguments};

    { // write header line
        std::string line{"#filenames"};
        for (auto const & filenames : index.bin_path())
        {
            line += '\t';
            for (auto const & filename : filenames)
            {
                line += filename;
                line += ';';
            }
        }
        line += '\n';
        synced_out.write(line);
    }

    std::vector<std::string> filenames{};

    auto worker = [&](size_t const start, size_t const end)
    {
        auto counter = index.ibf().template counting_agent<uint32_t>();

        std::string result_string{};

        for (auto && filename : filenames | seqan3::views::slice(start, end))
        {
            result_string.clear();
            result_string += filename;

            my_priority_queue<uint64_t> sketch{};
            for (auto && rec : seqan3::sequence_file_input<my_traits>{filename})
            {
                if (sketch.empty())
                    init_sketch(rec.sequence(), options.kmer_size, options.sketch_size, sketch);
                else
                    add_to_sketch(rec.sequence(), options.kmer_size, sketch);
            }

            auto & result = counter.bulk_count(sketch.get_underlying_container());

            for (size_t i = 0; i < result.size(); ++i)
            {
                auto const dist = compute_distance(result[i],
                                                   options.sketch_size,
                                                   options.fpr,
                                                   options.sizes.at(filename),
                                                   options.sizes.at(index.bin_path()[i][0]));

                result_string += '\t';
                result_string += std::to_string(dist);
            }

            result_string += '\n';
            synced_out.write(result_string);
        }
    };

    for (auto && chunked_files : options.files | seqan3::views::chunk((1ULL << 20) * 10))
    {
        filenames.clear();
        std::ranges::move(chunked_files, std::back_inserter(filenames));

        raptor::do_parallel(worker, filenames.size(), options.threads);
    }
}
