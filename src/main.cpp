#include <sstream>

#include <seqan3/argument_parser/all.hpp>

#include "search.hpp"
#include "jaqquard_dist.hpp"

int parse_command_line(smash_options & options, int const argc, char const * const * argv)
{
    seqan3::argument_parser parser{"smash", argc, argv};

    // Parser
    parser.info.author = "SeqAn-Team"; // give parser some infos
    parser.info.version = "1.0.0";
    parser.add_option(options.input_file, 'i', "input", "Please provide a file with one line one file each.");
    parser.add_option(options.index_file, 'x', "index", "Please provide an index file.");
    parser.add_option(options.output_file, 'o', "output", "The file for the distances matrix");
    parser.add_option(options.kmer_size, 'k', "kemr-size", "The kmer size.");
    parser.add_option(options.sketch_size, 's', "sketch-size", "The sketch size.");
    parser.add_option(options.threads, 't', "threads", "The number of threads to use.");
    parser.add_flag(options.no_sketching, 'd', "disable-sketching", "this will compute the true jaqquard distance.");

    try
    {
        parser.parse();                                                  // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        std::cerr << "Parsing error. " << ext.what() << "\n"; // give error message
        return -1;
    }

    return 0;
}

int main(int argc, char ** argv)
{
    smash_options options{};
    parse_command_line(options, argc, argv);

    std::string line;
    std::ifstream in{options.input_file};

    while (std::getline(in, line))
        options.files.push_back(line);

    if (options.no_sketching)
        jaqquard_dist(options);
    else
        search(options);

    return 0;
}
