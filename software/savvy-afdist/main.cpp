#include <savvy/reader.hpp>
#include <vector>
#include <thread>
#include <iostream>
#include <string>
#include <fstream>

void process_slice(const std::string& filename, uint64_t start, uint64_t end, std::vector<double>& thread_af, std::vector<int>& thread_hets, std::vector<int>& thread_homs) {
    savvy::reader f(filename);
    f.reset_bounds({start, end});
    savvy::variant var;
    std::vector<int> geno;

    while (f >> var) {
        var.get_format("GT", geno);
        int total_alleles = 0;
        int alt_alleles = 0;
        int het_count = 0;
        int hom_count = 0;

        for (size_t i = 0; i < geno.size(); i += 2) {
            int allele1 = geno[i];
            int allele2 = geno[i + 1];
            total_alleles += 2;
            alt_alleles += allele1 + allele2;

            if (allele1 != allele2) {
                het_count++;
            } else if (allele1 == 1 && allele2 == 1) {
                hom_count++;
            }
        }

        double allele_freq = total_alleles > 0 ? static_cast<double>(alt_alleles) / total_alleles : 0.0;
        thread_af.push_back(allele_freq);
        thread_hets.push_back(het_count);
        thread_homs.push_back(hom_count);
    }
}

// Use a lambda expression for finding the bin index
auto findBinIndex = [](const std::vector<double>& bins, double value) -> int {
    if (value < bins.front()) return 0;
    if (value >= bins.back()) return static_cast<int>(bins.size()) - 2;

    auto it = std::lower_bound(bins.begin(), bins.end(), value);
    return static_cast<int>((it == bins.end() || *it != value) ? it - bins.begin() - 1 : it - bins.begin());
};

// Function to check if a file exists
bool fileExists(const std::string& filename) {
    std::ifstream file(filename);
    return file.good();
}

int main(int argc, char* argv[]) {
    int num_threads = std::thread::hardware_concurrency(); // Default to the number of cores
    std::string filename;
    int total_records = -1; // Initialize total_records to an invalid state

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--threads") {
            if (i + 1 < argc) {
                num_threads = std::stoi(argv[++i]);
            }
        } else if (arg == "--sav-file") {
            if (i + 1 < argc) {
                filename = argv[++i];
            }
        } else if (arg == "--num-variants") {
            if (i + 1 < argc) {
                total_records = std::stoi(argv[++i]);
            }
        }
    }

    std::cout << "Number of threads: " << num_threads << std::endl;
    
    if (filename.empty()) {
        std::cerr << "Error: --sav-file argument is required" << std::endl;
        return 1;
    }

    if (!fileExists(filename)) {
        std::cerr << "Error: File '" << filename << "' does not exist" << std::endl;
        return 1;
    }

    if (total_records <= 0) {
        std::cerr << "Error: --num-variants argument is required and must be positive" << std::endl;
        return 1;
    }

    std::vector<std::thread> threads;
    std::vector<std::vector<double>> all_af(num_threads);
    std::vector<std::vector<int>> all_hets(num_threads);
    std::vector<std::vector<int>> all_homs(num_threads);

    // Divide the total number of records equally among the threads
    int records_per_thread = total_records / num_threads;

    for (int i = 0; i < num_threads; ++i) {
        int start = i * records_per_thread;
        int end = (i + 1) * records_per_thread;
        threads.emplace_back(process_slice, filename, start, end, std::ref(all_af[i]), std::ref(all_hets[i]), std::ref(all_homs[i]));
    }

    for (auto& t : threads) {
        t.join();
    }

    // Merge results
    std::vector<double> af;
    std::vector<int> hets;
    std::vector<int> homs;
    for (int i = 0; i < num_threads; ++i) {
        af.insert(af.end(), all_af[i].begin(), all_af[i].end());
        hets.insert(hets.end(), all_hets[i].begin(), all_hets[i].end());
        homs.insert(homs.end(), all_homs[i].begin(), all_homs[i].end());
    }

    std::cout << "Number of variants: " << af.size() << std::endl;

    const int num_bins = 10;
    std::vector<double> bins(num_bins + 1);
    for (int i = 0; i <= num_bins; ++i) {
        bins[i] = i * (1.0 / num_bins);
    }
    bins.back() += 0.01; // Adjust the last bin

    std::vector<int> het_bins(num_bins, 0), hom_bins(num_bins, 0), total_counts(num_bins, 0);

    int index = 0; 
    for (auto& freq : af) {
        double pRA = 2 * freq * (1 - freq);
        double pAA = freq * freq;

        int binIndex = findBinIndex(bins, pRA);
        het_bins[binIndex] += hets[index]; 

        binIndex = findBinIndex(bins, pAA);
        hom_bins[binIndex] += homs[index]; 

        index++; 
    }

    // Summing het and hom counts
    std::transform(het_bins.begin(), het_bins.end(), hom_bins.begin(), total_counts.begin(), std::plus<int>());

    // Output the results
    for (int i = 0; i < num_bins; ++i) {
        std::cout << "Bin " << i << " (" << bins[i] << " - " << bins[i + 1] << "): " << total_counts[i] << std::endl;
    }

    return 0;
}