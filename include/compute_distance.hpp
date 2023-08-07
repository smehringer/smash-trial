#pragma once

// computes the Jaqquard Index Value (sorry for the name)
// A = what I sketch/search from
// B = whats stored in my Bloom filter
double compute_distance(uint64_t const count,
                        uint64_t const sketch_size,
                        double const fpr,
                        uint64_t const size_of_A,
                        uint64_t const size_of_B)
{
    // Todo: is the -fpr always correct? isn't the impact stronger when JI is low ?

    // Containement estimate
    // Paper: C_est = ( Y^k / k ) − p
    // Here : C_est = ( count / sketch_size ) − fpr
    double const C_est = (static_cast<double>(count) / static_cast<double>(sketch_size)) - fpr;

    // J_est = |A|C_est / ( |A|+|B|−|A|C_est )
    double const A_C_est = size_of_A * C_est;
    double const J_est = A_C_est / (size_of_A + size_of_B - A_C_est);

    return J_est;
}

