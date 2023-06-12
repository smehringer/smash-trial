#include <charconv>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <ranges>

#include <robin_hood.h>

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

int main(int /* argc */, char ** argv)
{
    std::ifstream fin{argv[1]};

    std::vector<std::string> ids = read_column_names(fin);
    std::vector<size_t> const col_permutation = get_permutation(ids);

    std::vector<std::string> row_names;
    std::vector<std::string> rows;

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
        std::vector<std::string> columns;
        while (it != splitted_line.end())
        {
            number.clear();
            std::ranges::copy(*it, std::back_inserter(number));
            columns.push_back(number);
            ++it;
        }

        if (col_permutation.size() != columns.size())
            throw std::runtime_error{"permutation size:" + std::to_string(col_permutation.size()) +
                                     "columns size: " + std::to_string(columns.size())};

        apply_permutation(col_permutation, columns);

        // put together row
        name.reserve(line.size());
        for (auto const & col : columns)
        {
            name.push_back('\t');
            name += col;
        }
        name.push_back('\n');

        rows.push_back(name);
    }

    apply_permutation(col_permutation, ids);
    for (auto const & id : ids) // write header line
        std::cout << '\t' << id;
    std::cout << '\n';

    auto const & row_permutation = get_permutation(row_names);
    apply_permutation(row_permutation, rows);
    for (auto const & row : rows)
        std::cout << row;
}
