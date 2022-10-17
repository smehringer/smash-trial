#include <seqan3/search/views/kmer_hash.hpp>

#include <raptor/adjust_seed.hpp>
#include <raptor/dna4_traits.hpp>
#include <raptor/search/do_parallel.hpp>
#include <raptor/search/load_index.hpp>
#include <raptor/search/sync_out.hpp>

#include "jaqquard_dist.hpp"
#include "options.hpp"

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4; // instead of dna5
};

void jaqquard_dist(smash_options const & options)
{
    auto index = raptor::raptor_index<raptor::index_structure::hibf>{};

    auto get_kmers = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{options.kmer_size}});

    double index_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

    auto cereal_worker = [&]()
    {
        raptor::load_index(index, options.index_file, index_io_time);
    };
    auto cereal_handle = std::async(std::launch::async, cereal_worker);

    std::vector<std::string> filenames_chunk_buffer{};
    std::vector<uint64_t> sizes{};

    raptor::sync_out synced_out{options.output_file};

    cereal_handle.wait();

    std::cerr << "Writing header line and computing reference sizes..." << std::endl;
    { // write header line
        std::string line{"#filenames"};
        for (auto const & filenames : index.bin_path())
        {
            line += '\t';
            std::set<uint64_t> hashes{};
            for (auto const & filename : filenames)
            {
                for (auto && rec : seqan3::sequence_file_input<my_traits>{filename})
                    for (auto const hash : rec.sequence() | get_kmers)
                        hashes.insert(hash);

                line += filename;
                line += ';';
            }
            sizes.push_back(hashes.size());
        }
        line += '\n';
        synced_out << line;
    }

    auto worker = [&](size_t const start, size_t const end)
    {
        auto counter = index.ibf().template counting_agent<uint32_t>();

        std::string result_string{};

        for (auto && filename : filenames_chunk_buffer | seqan3::views::slice(start, end))
        {
            result_string.clear();
            result_string += filename;

            std::set<uint64_t> hashes{};

            for (auto && rec : seqan3::sequence_file_input<my_traits>{filename})
                for (auto const hash : rec.sequence() | get_kmers)
                    hashes.insert(hash);

            // For all hashes computed for current `filename` count their occurence for each user bin in the HIBF
            auto & result = counter.bulk_count(hashes);

            for (size_t i = 0; i < result.size(); ++i)
            {
                /* A intersect B    =                  result[i]               // #shared-hashes
                 * -------------    =   ------------------------------------
                 *   A union B      =   sizes[i] + hashes.size() - result[i]  // #hashes-A + #hashes-B - #shared-hashes
                 */
                auto const dist = static_cast<double>(result[i]) / (sizes[i] + hashes.size() - result[i]);

                result_string += '\t';
                result_string += std::to_string(dist);
            }

            result_string += '\n';
            synced_out.write(result_string);
        }
    };

    std::cerr << "Computing distances..." << std::endl;
    for (auto && chunk : options.files | seqan3::views::chunk((1ULL << 20) * 10))
    {
        filenames_chunk_buffer.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunk, std::back_inserter(filenames_chunk_buffer));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

        raptor::do_parallel(worker, filenames_chunk_buffer.size(), options.threads, compute_time);
    }

    // GCOVR_EXCL_START
    if (options.write_time)
    {
        std::filesystem::path file_path{options.output_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "Index I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed << std::setprecision(2) << index_io_time << '\t' << reads_io_time << '\t'
                    << compute_time;
    }
    // GCOVR_EXCL_STOP
}
