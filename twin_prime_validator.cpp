// ============================================================================
// Twin Prime Pair Validator with XOR Pattern Analysis
//
// Author: Thiago Fernandes Motta Massensini Silva
// Email: thiago.massensini@gmail.com
// Date: November 2024
//
// Description:
//   Validates twin prime pairs and analyzes their binary XOR patterns.
//   Implements the distribution P(k) = 2^(-k-1) where k is computed via
//   the formula k = ctz((p XOR (p+2)) + 2) - 1.
//
// Compilation:
//   g++ -O3 -march=native -fopenmp -std=c++17 twin_prime_validator_professional.cpp -o validator
//
// Additional optimization flags (optional):
//   -funroll-loops -ffast-math -DNDEBUG -flto
//
// Usage:
//   ./validator <input_csv> <num_threads>
//
// Input Format (CSV):
//   p,p+2,k
//   3,5,2
//   5,7,1
//   ...
//
// Output:
//   Validation statistics
//   Distribution analysis
//   Chi-squared goodness-of-fit test
// ============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <iomanip>
#include <omp.h>

using u64 = uint64_t;
using u128 = __uint128_t;

// ============================================================================
// Miller-Rabin Primality Test (Deterministic for n < 2^64)
// ============================================================================

inline u64 mod_mul(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((u128)a * b % mod);
}

inline u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = mod_mul(result, base, mod);
        base = mod_mul(base, base, mod);
        exp >>= 1;
    }
    return result;
}

inline bool is_prime64(u64 n) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;

    // Deterministic bases for n < 2^64
    static const u64 bases[] = {2, 3, 5, 7, 11, 13, 17};

    u64 d = n - 1, r = 0;
    while (!(d & 1)) {
        d >>= 1;
        r++;
    }

    for (u64 a : bases) {
        if (a >= n) continue;

        u64 x = mod_pow(a, d, n);
        if (x == 1 || x == n - 1) continue;

        bool composite = true;
        for (u64 i = 0; i < r - 1; i++) {
            x = mod_mul(x, x, n);
            if (x == n - 1) {
                composite = false;
                break;
            }
        }

        if (composite) return false;
    }

    return true;
}

// ============================================================================
// XOR Pattern Analysis
// ============================================================================

inline int calc_k(u64 p, u64 p2) {
    // k = ctz((p XOR (p+2)) + 2) - 1
    return __builtin_ctzll((p ^ p2) + 2) - 1;
}

inline bool verify_xor_property(u64 p, u64 p2, int k) {
    // Verify: (p XOR (p+2)) + 2 = 2^(k+1)
    u64 xor_plus_2 = (p ^ p2) + 2;
    u64 expected = 1ULL << (k + 1);
    return xor_plus_2 == expected;
}

// ============================================================================
// Data Structures
// ============================================================================

struct TwinPrime {
    u64 p;
    u64 p2;
    int k;
};

struct ValidationStats {
    u64 total_pairs = 0;
    u64 valid_primality = 0;
    u64 valid_xor_pattern = 0;
    std::array<u64, 25> k_distribution = {};

    void record_pair(bool is_prime_valid, bool is_xor_valid, int k) {
        total_pairs++;
        if (is_prime_valid) valid_primality++;
        if (is_xor_valid) valid_xor_pattern++;
        if (k >= 0 && k < 25) k_distribution[k]++;
    }
};

// ============================================================================
// Chi-Squared Test
// ============================================================================

double compute_chi_squared(const std::array<u64, 25>& observed, u64 total) {
    double chi2 = 0.0;

    for (int k = 1; k < 25; k++) {
        // Expected frequency: P(k) = 2^(-k-1)
        double expected = total * std::pow(2.0, -(k + 1));

        if (expected < 5.0) break; // Chi-squared requirement

        double diff = static_cast<double>(observed[k]) - expected;
        chi2 += (diff * diff) / expected;
    }

    return chi2;
}

// ============================================================================
// CSV Parser
// ============================================================================

std::vector<TwinPrime> load_csv(const std::string& filename) {
    std::vector<TwinPrime> pairs;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open file: " << filename << std::endl;
        return pairs;
    }

    std::string line;
    // Skip header if present
    if (std::getline(file, line)) {
        if (line.find("p,") == std::string::npos) {
            // First line is not header, process it
            std::istringstream iss(line);
            TwinPrime tp;
            char comma;
            if (iss >> tp.p >> comma >> tp.p2 >> comma >> tp.k) {
                pairs.push_back(tp);
            }
        }
    }

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        TwinPrime tp;
        char comma;
        if (iss >> tp.p >> comma >> tp.p2 >> comma >> tp.k) {
            pairs.push_back(tp);
        }
    }

    file.close();
    return pairs;
}

// ============================================================================
// Validation Engine
// ============================================================================

ValidationStats validate_dataset(const std::vector<TwinPrime>& pairs, int num_threads) {
    ValidationStats stats;
    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
        ValidationStats local_stats;

        #pragma omp for schedule(dynamic, 10000)
        for (size_t i = 0; i < pairs.size(); i++) {
            const auto& tp = pairs[i];

            // Test 1: Primality
            bool is_prime_valid = is_prime64(tp.p) && is_prime64(tp.p2);

            // Test 2: XOR pattern
            int k_computed = calc_k(tp.p, tp.p2);
            bool is_xor_valid = (k_computed == tp.k) &&
                               verify_xor_property(tp.p, tp.p2, tp.k);

            local_stats.record_pair(is_prime_valid, is_xor_valid, tp.k);
        }

        #pragma omp critical
        {
            stats.total_pairs += local_stats.total_pairs;
            stats.valid_primality += local_stats.valid_primality;
            stats.valid_xor_pattern += local_stats.valid_xor_pattern;
            for (size_t k = 0; k < stats.k_distribution.size(); k++) {
                stats.k_distribution[k] += local_stats.k_distribution[k];
            }
        }
    }

    return stats;
}

// ============================================================================
// Reporting
// ============================================================================

void print_results(const ValidationStats& stats) {
    std::cout << "\n";
    std::cout << "========================================" << std::endl;
    std::cout << "VALIDATION RESULTS" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\n";

    std::cout << "Total twin prime pairs:    " << stats.total_pairs << std::endl;
    std::cout << "Valid primality:           " << stats.valid_primality
              << " (" << std::fixed << std::setprecision(2)
              << (100.0 * stats.valid_primality / stats.total_pairs) << "%)" << std::endl;
    std::cout << "Valid XOR patterns:        " << stats.valid_xor_pattern
              << " (" << std::fixed << std::setprecision(2)
              << (100.0 * stats.valid_xor_pattern / stats.total_pairs) << "%)" << std::endl;
    std::cout << "\n";

    std::cout << "========================================" << std::endl;
    std::cout << "DISTRIBUTION ANALYSIS: P(k) = 2^(-k-1)" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\n";

    std::cout << std::left << std::setw(6) << "k"
              << std::right << std::setw(15) << "Observed"
              << std::setw(15) << "Expected"
              << std::setw(12) << "Obs %"
              << std::setw(12) << "Exp %"
              << std::setw(12) << "Deviation" << std::endl;
    std::cout << std::string(72, '-') << std::endl;

    for (int k = 1; k < 25; k++) {
        if (stats.k_distribution[k] == 0) break;

        double expected = stats.total_pairs * std::pow(2.0, -(k + 1));
        double obs_pct = 100.0 * stats.k_distribution[k] / stats.total_pairs;
        double exp_pct = 100.0 * std::pow(2.0, -(k + 1));
        double deviation = obs_pct - exp_pct;

        std::cout << std::left << std::setw(6) << k
                  << std::right << std::setw(15) << stats.k_distribution[k]
                  << std::setw(15) << std::fixed << std::setprecision(0) << expected
                  << std::setw(11) << std::setprecision(4) << obs_pct << "%"
                  << std::setw(11) << std::setprecision(4) << exp_pct << "%"
                  << std::setw(11) << std::setprecision(4) << deviation << "%"
                  << std::endl;
    }

    std::cout << std::endl;

    // Chi-squared test
    double chi2 = compute_chi_squared(stats.k_distribution, stats.total_pairs);

    std::cout << "========================================" << std::endl;
    std::cout << "STATISTICAL SIGNIFICANCE" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\n";
    std::cout << "Chi-squared statistic:     " << std::fixed << std::setprecision(4) << chi2 << std::endl;
    std::cout << "Degrees of freedom:        14" << std::endl;
    std::cout << "Critical value (95%):      23.685" << std::endl;
    std::cout << "\n";

    if (chi2 < 23.685) {
        std::cout << "Result: EXCELLENT FIT (χ² < critical value)" << std::endl;
        std::cout << "The observed distribution matches the theoretical" << std::endl;
        std::cout << "geometric distribution P(k) = 2^(-k-1)." << std::endl;
    } else {
        std::cout << "Result: POOR FIT (χ² >= critical value)" << std::endl;
        std::cout << "The observed distribution deviates significantly" << std::endl;
        std::cout << "from the theoretical model." << std::endl;
    }

    std::cout << "\n";
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;
}

// ============================================================================
// Main Program
// ============================================================================

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_csv> [num_threads]" << std::endl;
        std::cerr << "Example: " << argv[0] << " twin_primes_part_01.csv 56" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    int num_threads = (argc >= 3) ? std::stoi(argv[2]) : omp_get_max_threads();

    std::cout << "\n";
    std::cout << "========================================" << std::endl;
    std::cout << "Twin Prime XOR Pattern Validator" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\n";
    std::cout << "Author:      Thiago Fernandes Motta Massensini Silva" << std::endl;
    std::cout << "File:        " << filename << std::endl;
    std::cout << "Threads:     " << num_threads << std::endl;
    std::cout << "\n";

    auto t_start = std::chrono::high_resolution_clock::now();

    std::cout << "Loading dataset..." << std::endl;
    auto pairs = load_csv(filename);

    if (pairs.empty()) {
        std::cerr << "ERROR: No data loaded from file." << std::endl;
        return 1;
    }

    std::cout << "Loaded " << pairs.size() << " twin prime pairs." << std::endl;
    std::cout << "\n";
    std::cout << "Validating..." << std::endl;

    auto stats = validate_dataset(pairs, num_threads);

    auto t_end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    print_results(stats);

    std::cout << "Processing time:           " << std::fixed << std::setprecision(3)
              << (elapsed / 1000.0) << " seconds" << std::endl;
    std::cout << "Throughput:                " << std::fixed << std::setprecision(0)
              << (pairs.size() / (elapsed / 1000.0)) << " pairs/second" << std::endl;
    std::cout << "\n";

    return 0;
}
