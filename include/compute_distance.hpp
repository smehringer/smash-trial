#pragma once

double compute_distance(uint64_t const count, uint64_t const sketch_size)
{
    return static_cast<double>(count) / sketch_size;
}
