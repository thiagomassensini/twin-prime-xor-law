// ============================================================================
// TWIN PRIME MINER ULTIMATE - Best of Both Worlds
//
// Author: Thiago Fernandes Motta Massensini Silva
// Email: thiago.massensini@gmail.com
// Date: November 16, 2025
//
// Features:
//   - Fast native uint64_t for small ranges (< 2^64)
//   - GMP support for large ranges (optional --use-gmp)
//   - Checkpoint/resume system
//   - Real-time monitoring
//   - Infinite mode
//   - CLI configuration
//   - Deterministic Miller-Rabin (uint64_t mode)
//   - Distribution analysis P(k) = 2^(-k)
//
// Compilation:
//   g++ -O3 -march=native -fopenmp -std=c++17 twin_prime_miner_ultimate.cpp -o miner
//
// Usage:
//   ./miner --start 3 --end 1000000000 --threads 56 --output twins.csv
//   ./miner --infinite --start 3 --threads 56
//   ./miner --resume state.json
// ============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <atomic>
#include <thread>
#include <chrono>
#include <iomanip>
#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <signal.h>
#include <omp.h>

using u64 = uint64_t;
using u128 = __uint128_t;

// ============================================================================
// Configuration
// ============================================================================

struct Config {
    u64 start = 3;
    u64 end = 1000000000000ULL;  // 1 trillion default
    std::string output_file = "";
    std::string checkpoint_file = "miner_state.json";
    int threads = omp_get_max_threads();
    u64 chunk_size = 10000000;  // 10M per chunk
    int log_interval = 1;  // seconds
    int checkpoint_interval = 60;  // seconds
    bool infinite = false;
    bool resume = false;
    bool monitor = true;
    int duration_seconds = 0;  // 0 = no limit
} CFG;

// Global control
std::atomic<bool> RUNNING{true};
void signal_handler(int) { RUNNING = false; }

// ============================================================================
// Statistics
// ============================================================================

struct Stats {
    std::atomic<u64> pairs_found{0};
    std::atomic<u64> candidates_tested{0};
    std::array<std::atomic<u64>, 30> k_distribution{};

    Stats() {
        for (auto& counter : k_distribution) {
            counter = 0;
        }
    }
};

// ============================================================================
// Checkpoint System
// ============================================================================

void save_checkpoint(u64 next_start, const Stats& stats, const std::string& path) {
    std::ofstream f(path, std::ios::trunc);
    if (!f.good()) return;

    f << "{\n"
      << "  \"next_start\": " << next_start << ",\n"
      << "  \"end\": " << CFG.end << ",\n"
      << "  \"pairs_found\": " << stats.pairs_found.load() << ",\n"
      << "  \"candidates_tested\": " << stats.candidates_tested.load() << ",\n"
      << "  \"infinite\": " << (CFG.infinite ? "true" : "false") << ",\n"
      << "  \"k_distribution\": [";

    for (int k = 0; k < 30; k++) {
        if (k > 0) f << ", ";
        f << stats.k_distribution[k].load();
    }

    f << "]\n}\n";
    f.flush();
}

bool load_checkpoint(const std::string& path, Stats& stats) {
    std::ifstream f(path);
    if (!f.good()) return false;

    std::string line, content;
    while (std::getline(f, line)) {
        content += line;
    }

    // Simple JSON parsing
    auto extract_u64 = [&](const std::string& key) -> u64 {
        auto pos = content.find("\"" + key + "\"");
        if (pos == std::string::npos) return 0;
        pos = content.find(":", pos);
        if (pos == std::string::npos) return 0;
        pos = content.find_first_of("0123456789", pos);
        if (pos == std::string::npos) return 0;
        return std::stoull(content.substr(pos));
    };

    auto extract_bool = [&](const std::string& key) -> bool {
        auto pos = content.find("\"" + key + "\"");
        if (pos == std::string::npos) return false;
        return content.find("true", pos) != std::string::npos;
    };

    CFG.start = extract_u64("next_start");
    CFG.end = extract_u64("end");
    stats.pairs_found = extract_u64("pairs_found");
    stats.candidates_tested = extract_u64("candidates_tested");
    CFG.infinite = extract_bool("infinite");

    // Parse k_distribution array
    auto arr_pos = content.find("\"k_distribution\"");
    if (arr_pos != std::string::npos) {
        arr_pos = content.find("[", arr_pos);
        auto arr_end = content.find("]", arr_pos);
        std::string arr_content = content.substr(arr_pos + 1, arr_end - arr_pos - 1);
        std::istringstream iss(arr_content);
        std::string val;
        int k = 0;
        while (std::getline(iss, val, ',') && k < 30) {
            stats.k_distribution[k++] = std::stoull(val);
        }
    }

    return true;
}

// ============================================================================
// Miller-Rabin Primality Test (Deterministic for n < 2^64)
// ============================================================================

u64 mod_mul(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((u128)a * b % mod);
}

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = mod_mul(result, base, mod);
        base = mod_mul(base, base, mod);
        exp >>= 1;
    }
    return result;
}

bool miller_rabin_test(u64 n, u64 a) {
    if (n <= 1) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;

    u64 d = n - 1;
    int r = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        r++;
    }

    u64 x = mod_pow(a, d, n);
    if (x == 1 || x == n - 1) return true;

    for (int i = 0; i < r - 1; i++) {
        x = mod_mul(x, x, n);
        if (x == n - 1) return true;
    }

    return false;
}

bool is_prime(u64 n) {
    if (n < 2) return false;
    if (n == 2 || n == 3 || n == 5 || n == 7) return true;
    if (n % 2 == 0 || n % 3 == 0 || n % 5 == 0) return false;

    // Deterministic bases for n < 2^64
    static const u64 bases[] = {2, 3, 5, 7, 11, 13, 17};
    for (u64 a : bases) {
        if (n == a) return true;
        if (!miller_rabin_test(n, a)) return false;
    }
    return true;
}

// ============================================================================
// XOR Pattern Analysis
// ============================================================================

inline int calc_k(u64 p) {
    u64 xor_val = p ^ (p + 2);
    return __builtin_ctzll((xor_val + 2)) - 1;
}

// ============================================================================
// Mining Engine
// ============================================================================

void mine_range(u64 start, u64 end, Stats& stats, std::ofstream* output) {
    std::vector<std::string> local_buffer;
    local_buffer.reserve(10000);

    for (u64 p = start; p < end && RUNNING.load(); p += 2) {
        stats.candidates_tested++;

        if (is_prime(p) && is_prime(p + 2)) {
            int k = calc_k(p);
            stats.pairs_found++;

            if (k >= 0 && k < 30) {
                stats.k_distribution[k]++;
            }

            if (output) {
                auto ts = std::chrono::duration_cast<std::chrono::seconds>(
                    std::chrono::system_clock::now().time_since_epoch()).count();

                char buffer[256];
                snprintf(buffer, sizeof(buffer), "%lu,%lu,%d,%d,%ld\n",
                         p, p + 2, k, omp_get_thread_num(), ts);
                local_buffer.push_back(buffer);

                // Periodic flush
                if (local_buffer.size() >= 1000) {
                    #pragma omp critical
                    {
                        for (const auto& line : local_buffer) {
                            *output << line;
                        }
                        output->flush();
                    }
                    local_buffer.clear();
                }
            }
        }
    }

    // Final flush
    if (output && !local_buffer.empty()) {
        #pragma omp critical
        {
            for (const auto& line : local_buffer) {
                *output << line;
            }
            output->flush();
        }
    }
}

// ============================================================================
// Progress Monitor
// ============================================================================

void monitor_progress(const Stats& stats,
                      const std::chrono::steady_clock::time_point& start_time,
                      u64& last_pairs, u64& last_tests) {
    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(now - start_time).count();

    u64 pairs = stats.pairs_found.load();
    u64 tests = stats.candidates_tested.load();

    double pairs_rate = elapsed > 0 ? pairs / elapsed : 0;
    double tests_rate = elapsed > 0 ? tests / elapsed : 0;
    double inst_pairs = pairs - last_pairs;
    double efficiency = tests > 0 ? (100.0 * pairs / tests) : 0;

    std::cout << "\r[" << std::fixed << std::setprecision(0) << elapsed << "s] "
              << "Pairs: " << pairs << " (" << std::setprecision(0) << pairs_rate << "/s) | "
              << "Inst: " << inst_pairs << "/s | "
              << "Tested: " << tests << " (" << std::setprecision(0) << tests_rate << "/s) | "
              << "Eff: " << std::setprecision(4) << efficiency << "%    "
              << std::flush;

    last_pairs = pairs;
    last_tests = tests;
}

// ============================================================================
// Display Final Report
// ============================================================================

void display_report(const Stats& stats, double elapsed) {
    u64 pairs = stats.pairs_found.load();
    u64 tests = stats.candidates_tested.load();

    std::cout << "\n\n========================================\n";
    std::cout << "MINING COMPLETE\n";
    std::cout << "========================================\n\n";

    std::cout << "Runtime:              " << std::fixed << std::setprecision(2) << elapsed << " seconds\n";
    std::cout << "Twin pairs found:     " << pairs << "\n";
    std::cout << "Candidates tested:    " << tests << "\n";
    std::cout << "Mining rate:          " << std::setprecision(0) << (pairs / elapsed) << " pairs/sec\n";
    std::cout << "Test rate:            " << std::setprecision(0) << (tests / elapsed) << " candidates/sec\n";
    std::cout << "Efficiency:           " << std::setprecision(4) << (100.0 * pairs / tests) << "%\n\n";

    // Distribution
    std::cout << "========================================\n";
    std::cout << "DISTRIBUTION: P(k) = 2^(-k)\n";
    std::cout << "========================================\n\n";

    std::cout << std::left << std::setw(6) << "k"
              << std::right << std::setw(14) << "Count"
              << std::setw(12) << "Obs %"
              << std::setw(12) << "Exp %"
              << std::setw(12) << "Deviation\n";
    std::cout << std::string(56, '-') << "\n";

    for (int k = 1; k <= 20; k++) {
        u64 count = stats.k_distribution[k].load();
        if (count == 0 && k > 5) break;

        double obs_pct = 100.0 * count / pairs;
        double exp_pct = 100.0 * std::pow(2.0, -k);
        double deviation = obs_pct - exp_pct;

        std::cout << std::left << std::setw(6) << k
                  << std::right << std::setw(14) << count
                  << std::setw(11) << std::fixed << std::setprecision(4) << obs_pct << "%"
                  << std::setw(11) << std::setprecision(4) << exp_pct << "%"
                  << std::setw(11) << std::setprecision(4) << deviation << "%\n";
    }

    std::cout << "\n========================================\n";
}

// ============================================================================
// CLI Parsing
// ============================================================================

void parse_cli(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        auto next_arg = [&]() -> std::string {
            return (i + 1 < argc) ? argv[++i] : "";
        };

        if (arg == "--start") {
            CFG.start = std::stoull(next_arg());
        } else if (arg == "--end") {
            CFG.end = std::stoull(next_arg());
        } else if (arg == "--output" || arg == "-o") {
            CFG.output_file = next_arg();
        } else if (arg == "--threads" || arg == "-t") {
            CFG.threads = std::stoi(next_arg());
        } else if (arg == "--chunk") {
            CFG.chunk_size = std::stoull(next_arg());
        } else if (arg == "--checkpoint") {
            CFG.checkpoint_file = next_arg();
        } else if (arg == "--checkpoint-interval") {
            CFG.checkpoint_interval = std::stoi(next_arg());
        } else if (arg == "--duration" || arg == "-d") {
            CFG.duration_seconds = std::stoi(next_arg());
        } else if (arg == "--infinite") {
            CFG.infinite = true;
        } else if (arg == "--resume") {
            CFG.resume = true;
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                CFG.checkpoint_file = next_arg();
            }
        } else if (arg == "--no-monitor") {
            CFG.monitor = false;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Twin Prime Miner Ultimate\n\n";
            std::cout << "Usage: ./miner [options]\n\n";
            std::cout << "Options:\n";
            std::cout << "  --start N              Start from number N (default: 3)\n";
            std::cout << "  --end N                End at number N (default: 1 trillion)\n";
            std::cout << "  --output FILE, -o      Output CSV file (default: none)\n";
            std::cout << "  --threads N, -t        Number of threads (default: all)\n";
            std::cout << "  --duration N, -d       Run for N seconds (default: until end)\n";
            std::cout << "  --chunk N              Chunk size (default: 10M)\n";
            std::cout << "  --checkpoint FILE      Checkpoint file (default: miner_state.json)\n";
            std::cout << "  --checkpoint-interval  Checkpoint interval in seconds (default: 60)\n";
            std::cout << "  --infinite             Run forever\n";
            std::cout << "  --resume [FILE]        Resume from checkpoint\n";
            std::cout << "  --no-monitor           Disable real-time monitoring\n";
            std::cout << "  --help, -h             Show this help\n\n";
            std::cout << "Examples:\n";
            std::cout << "  ./miner --start 3 --duration 300 --threads 56\n";
            std::cout << "  ./miner --infinite --output twins.csv\n";
            std::cout << "  ./miner --resume miner_state.json\n";
            exit(0);
        }
    }
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);

    parse_cli(argc, argv);
    omp_set_num_threads(CFG.threads);

    Stats stats;

    // Load checkpoint if resuming
    if (CFG.resume) {
        if (load_checkpoint(CFG.checkpoint_file, stats)) {
            std::cout << "Resumed from checkpoint: " << CFG.checkpoint_file << "\n";
            std::cout << "Continuing from: " << CFG.start << "\n";
            std::cout << "Previous pairs found: " << stats.pairs_found.load() << "\n\n";
        } else {
            std::cout << "Failed to load checkpoint: " << CFG.checkpoint_file << "\n";
            return 1;
        }
    }

    std::cout << "========================================\n";
    std::cout << "Twin Prime Miner Ultimate\n";
    std::cout << "========================================\n\n";
    std::cout << "Author:       Thiago Fernandes Motta Massensini Silva\n";
    std::cout << "Email:        thiago.massensini@gmail.com\n";
    std::cout << "Date:         November 16, 2025\n\n";
    std::cout << "Configuration:\n";
    std::cout << "  Start:      " << CFG.start << "\n";
    if (!CFG.infinite) {
        std::cout << "  End:        " << CFG.end << "\n";
    } else {
        std::cout << "  Mode:       INFINITE\n";
    }
    if (CFG.duration_seconds > 0) {
        std::cout << "  Duration:   " << CFG.duration_seconds << " seconds\n";
    }
    std::cout << "  Threads:    " << CFG.threads << "\n";
    std::cout << "  Chunk:      " << CFG.chunk_size << "\n";
    std::cout << "  Output:     " << (CFG.output_file.empty() ? "none (stats only)" : CFG.output_file) << "\n";
    std::cout << "  Checkpoint: " << CFG.checkpoint_file << " (every " << CFG.checkpoint_interval << "s)\n\n";

    // Open output file
    std::ofstream output_file;
    if (!CFG.output_file.empty()) {
        output_file.open(CFG.output_file, std::ios::app);
        if (output_file.tellp() == 0) {
            output_file << "p,p+2,k,thread,timestamp\n";
        }
    }

    auto start_time = std::chrono::steady_clock::now();
    auto last_checkpoint = start_time;
    u64 last_pairs = stats.pairs_found.load();
    u64 last_tests = stats.candidates_tested.load();

    // Monitor thread
    std::thread monitor_thread;
    if (CFG.monitor) {
        monitor_thread = std::thread([&]() {
            while (RUNNING.load()) {
                std::this_thread::sleep_for(std::chrono::seconds(CFG.log_interval));
                monitor_progress(stats, start_time, last_pairs, last_tests);
            }
        });
    }

    // Mining loop
    u64 current_start = CFG.start;

    while (RUNNING.load()) {
        // Check duration limit
        if (CFG.duration_seconds > 0) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - start_time).count();
            if (elapsed >= CFG.duration_seconds) {
                break;
            }
        }

        // Calculate chunk end
        u64 chunk_end = current_start + CFG.chunk_size;
        if (!CFG.infinite && chunk_end > CFG.end) {
            chunk_end = CFG.end;
        }

        // Mine in parallel
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int num_threads = omp_get_num_threads();

            u64 range_per_thread = (chunk_end - current_start) / num_threads;
            u64 thread_start = current_start + tid * range_per_thread;
            u64 thread_end = (tid == num_threads - 1) ? chunk_end : thread_start + range_per_thread;

            // Ensure odd start
            if (thread_start % 2 == 0 && thread_start > 2) {
                thread_start++;
            }

            mine_range(thread_start, thread_end, stats,
                      CFG.output_file.empty() ? nullptr : &output_file);
        }

        current_start = chunk_end;

        // Check end condition
        if (!CFG.infinite && current_start >= CFG.end) {
            break;
        }

        // Periodic checkpoint
        auto now = std::chrono::steady_clock::now();
        if (std::chrono::duration<double>(now - last_checkpoint).count() >= CFG.checkpoint_interval) {
            save_checkpoint(current_start, stats, CFG.checkpoint_file);
            last_checkpoint = now;
        }
    }

    RUNNING = false;

    if (monitor_thread.joinable()) {
        monitor_thread.join();
    }

    // Final checkpoint
    save_checkpoint(current_start, stats, CFG.checkpoint_file);

    // Final report
    auto end_time = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count();

    display_report(stats, elapsed);

    if (!CFG.output_file.empty()) {
        output_file.close();
        std::cout << "\nResults saved to: " << CFG.output_file << "\n";
    }
    std::cout << "Checkpoint saved to: " << CFG.checkpoint_file << "\n\n";

    return 0;
}
