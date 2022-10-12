#include <charconv>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <ranges>

#include <robin_hood.h>

int main(int /* argc */, char ** argv)
{
    std::ifstream fin{argv[1]};

    robin_hood::unordered_map<std::string, uint64_t> names;

    std::vector<std::vector<double>> matrix;

    std::string line;
    uint64_t counter{0};
    while (std::getline(fin, line))
    {
        auto splitted_line = line | std::views::split('\t');
        auto it = splitted_line.begin();

        std::string first{};
        std::ranges::copy(*it, std::back_inserter(first));
        ++it;
        std::string second{};
        std::ranges::copy(*it, std::back_inserter(second));
        ++it;
        ++it;
        ++it;
        std::string numbers_str{};
        std::ranges::copy(*it, std::back_inserter(numbers_str));

        double nominator{};
        double denominator{};

        auto res = std::from_chars(numbers_str.data(), numbers_str.data() + numbers_str.size(), nominator);
        ++res.ptr;
        res = std::from_chars(res.ptr, numbers_str.data() + numbers_str.size(), denominator);

        double dist = nominator / denominator;

        uint64_t column_index;
        uint64_t row_index;

        // find row_idx
        auto row = names.find(first);
        if (row == names.end())
        {
            names[first] = counter;
            row_index = counter;
            for (auto & row : matrix)
                row.push_back(0.);
            ++counter;
            matrix.push_back(std::vector<double>(counter, 0.));
        }
        else
        {
            row_index = row->second;
        }

        // find column index
        auto col = names.find(second);
        if (col == names.end())
        {
            names[first] = counter;
            column_index = counter;
            for (auto & row : matrix)
                row.push_back(0.);
            ++counter;
            matrix.push_back(std::vector<double>(counter, 0.));
        }
        else
        {
            column_index = col->second;
        }

        assert(matrix.size() > row_index);
        assert(matrix[row_index].size() > column_index);
        matrix[row_index][column_index] = dist;
    }

    std::vector<std::string> ids{};
    ids.reserve(names.size());
    for (auto & [name, idx] : names)
        ids.push_back(name);

    std::sort(ids.begin(), ids.end());

    for (auto const & id1 : ids) // write header line
        std::cout << '\t' << id1;
    std::cout << '\n';

    for (auto const & id1 : ids)
    {
        std::cout << id1;
        for (auto const & id2 : ids)
        {
            std::cout << '\t' << matrix[names.at(id1)][names.at(id2)];
        }
        std::cout << '\n';
    }
}
