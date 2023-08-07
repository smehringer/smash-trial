#include <charconv>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <ranges>

#include <robin_hood.h>

struct mash_data
{
    std::string first{};
    std::string second{};
    double dist{0.0};
};

mash_data parse_mash_line(std::string const & line)
{
    mash_data result{};
    auto splitted_line = line | std::views::split('\t');
    auto it = splitted_line.begin();

    std::ranges::copy(*it, std::back_inserter(result.first));
    ++it;
    std::ranges::copy(*it, std::back_inserter(result.second));
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

    result.dist = nominator / denominator;

    return result;
}

std::vector<mash_data> read_and_sort_native_mash_file(std::string const & filename)
{
    std::vector<mash_data> mash_lines;

    std::ifstream fin{filename};
    std::string line;

    while (std::getline(fin, line))
        mash_lines.push_back(parse_mash_line(line));

    // sort by second column and then by first column
    std::sort(mash_lines.begin(), mash_lines.end(),
              [](auto const & c1, auto const & c2)
              {
                  if (c1.second == c2.second)
                      return c1.first < c2.first;
                  else
                      return c1.second < c2.second;
              });

    return mash_lines;
}

void write_header_line(std::vector<mash_data> const & mash_lines)
{
    std::string const current = mash_lines[0].second;

    std::cout << "#filenames";

    for (mash_data const & mash_output : mash_lines)
    {
        if (current == mash_output.second)
            std::cout << '\t' << mash_output.first;
        else
            break;
    }
    std::cout << '\n';
}

int main(int /* argc */, char ** argv)
{
    auto const & mash_lines = read_and_sort_native_mash_file(argv[1]);
    assert(!mash_lines.empty());

    write_header_line(mash_lines);

    std::string current = mash_lines[0].second;
    std::cout << current;

    for (mash_data const & mash_output : mash_lines)
    {
        if (current == mash_output.second)
        {
            std::cout << '\t' << mash_output.dist;
        }
        else
        {
            std::cout << '\n' << mash_output.second << '\t' << mash_output.dist;
            current = mash_output.second;
        }
    }
    std::cout << '\n';
}
