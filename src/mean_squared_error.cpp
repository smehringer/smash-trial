#include <sstream>
#include <cmath>

#include <seqan3/argument_parser/all.hpp>

struct error_options
{
    std::string input_maxtrix_filename{};
    std::string truth_maxtrix_filename{};
};

template<typename random_access_range_type>
void apply_permutation(std::vector<size_t> const & permutation, random_access_range_type & permute_me)
{
    for (size_t i{0}; i < permutation.size(); ++i)
    {
        size_t swap_index = permutation[i];
        while (swap_index < i)
            swap_index = permutation[swap_index];

        std::swap(permute_me[i], permute_me[swap_index]);
    }
}

int parse_command_line(error_options & options, int const argc, char const * const * argv)
{
    seqan3::argument_parser parser{"smash", argc, argv};

    // Parser
    parser.info.author = "SeqAn-Team"; // give parser some infos
    parser.info.version = "1.0.0";
    parser.add_option(options.input_maxtrix_filename, '\0', "input", "Please provide a file with a matrix.");
    parser.add_option(options.truth_maxtrix_filename, '\0', "truth", "Please provide the truth matrix of the same size as input.");

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

std::vector<std::string> read_column_names(std::ifstream & fin)
{
    std::vector<std::string> ids{};
    std::string line;

    // read header line for column names
    {
        std::getline(fin, line);

        auto splitted_line = line | std::views::split('\t');
        auto it = splitted_line.begin();
        ++it; // skip `#filenames` in the beginning

        std::string name{};
        while (it != splitted_line.end())
        {
            name.clear();
            std::ranges::copy(*it, std::back_inserter(name));
            if (name.back() == ';')
                name.pop_back();
            ids.push_back(name);
            ++it;
        }
    }

    return ids;
}

std::vector<size_t> get_permutation(std::vector<std::string> const & ids)
{
    std::vector<size_t> col_permutation{};
    col_permutation.reserve(ids.size());
    for (size_t i = 0; i < ids.size(); ++i)
        col_permutation.push_back(i);
    std::sort(col_permutation.begin(), col_permutation.end(), [&ids](auto const & i1, auto const & i2){ return ids[i1] < ids[i2]; });
    return col_permutation;
}

std::pair<std::vector<std::string>, std::vector<std::vector<double>>>
read_matrix(std::ifstream & fin, std::vector<size_t> const & col_permutation)
{
    std::vector<std::vector<double>> matrix;
    std::vector<std::string> row_names;
    std::string line;

    while (std::getline(fin, line))
    {
        auto splitted_line = line | std::views::split('\t');
        auto it = splitted_line.begin();

        std::string name{};
        std::ranges::copy(*it, std::back_inserter(name));
        row_names.push_back(name);
        ++it;

        std::string number{};
        std::vector<double> columns;
        while (it != splitted_line.end())
        {
            number.clear();
            double tmp;
            std::ranges::copy(*it, std::back_inserter(number));
            std::from_chars(number.data(), number.data() + number.size(), tmp);
            columns.push_back(tmp);
            ++it;
        }

        if (col_permutation.size() != columns.size())
            throw std::runtime_error{"permutation size:" + std::to_string(col_permutation.size()) +
                                     "columns size: " + std::to_string(columns.size())};

        apply_permutation(col_permutation, columns);

        matrix.push_back(std::move(columns));
    }

    return std::make_pair(row_names, matrix);
}

int main(int argc, char ** argv)
{
    error_options options{};
    parse_command_line(options, argc, argv);

    std::ifstream input_maxtrix_file{options.input_maxtrix_filename};
    std::ifstream truth_maxtrix_file{options.truth_maxtrix_filename};

    std::vector<std::string> input_columnnames = read_column_names(input_maxtrix_file);
    std::vector<std::string> truth_columnnames = read_column_names(truth_maxtrix_file);

    if (input_columnnames.size() != truth_columnnames.size())
        throw std::runtime_error{"ERROR: input and truth do not have the same number of columns."};

    std::vector<size_t> const input_col_permutation = get_permutation(input_columnnames);
    std::vector<size_t> const truth_col_permutation = get_permutation(truth_columnnames);

    for (size_t i = 0; i < input_columnnames.size(); ++i)
    {
        if (input_columnnames[input_col_permutation[i]] != truth_columnnames[truth_col_permutation[i]])
            throw std::runtime_error{"ERROR: Sorted columns of input and truth do not match. At position " +
                                     std::to_string(i) + " with input = " + input_columnnames[input_col_permutation[i]] +
                                     " truth = " + truth_columnnames[truth_col_permutation[i]]};
    }

    auto && [input_rownames, input_matrix] = read_matrix(input_maxtrix_file, input_col_permutation);
    auto && [truth_rownames, truth_matrix] = read_matrix(truth_maxtrix_file, truth_col_permutation);

    if (input_rownames.size() != truth_rownames.size())
        throw std::runtime_error{"ERROR: input and truth do not have the same number of rows."};

    auto const & input_row_permutation = get_permutation(input_rownames);
    auto const & truth_row_permutation = get_permutation(truth_rownames);

    for (size_t i = 0; i < input_rownames.size(); ++i)
    {
        if (input_rownames[input_row_permutation[i]] != truth_rownames[truth_row_permutation[i]])
            throw std::runtime_error{"ERROR: Sorted columns of input and truth do not match. At position " +
                                     std::to_string(i) + " with input = " + input_rownames[input_row_permutation[i]] +
                                     " truth = " + truth_rownames[truth_row_permutation[i]]};
    }

    apply_permutation(input_row_permutation, input_matrix);
    apply_permutation(truth_row_permutation, truth_matrix);

    // finally the matrixes are read in and sorted and seem to have the same rownames and column names
    // now we can compute the mean squared error:

    double sum{};
    for (size_t row_idx = 0; row_idx < input_rownames.size(); ++row_idx)
    {
        for (size_t column_idx = 0; column_idx < input_columnnames.size(); ++column_idx)
        {
            sum += std::pow(truth_matrix[row_idx][column_idx] - input_matrix[row_idx][column_idx], 2);
        }
    }

    std::cout << "SSE: " << sum << std::endl;

    return 0;
}
