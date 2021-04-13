#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include "io_utils.h"

using namespace std;

uint32_t kmer_count_len = 0;
map<uint32_t, uint32_t> kmer_inds;
uint32_t k_size;

u_int64_t rev_comp(u_int64_t x)
{
    u_int64_t res = x;

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;

    return (res >> (2 * (32 - k_size)));
}

void compute_kmer_inds()
{
    u_int64_t kmer_rc, ind = 0;

    for (u_int64_t kmer = 0; kmer < (u_int64_t)pow(4, k_size); kmer++)
    {
        kmer_rc = rev_comp(kmer);

        if (kmer_inds.find(kmer_rc) != kmer_inds.end())
        {
            kmer_inds[kmer] = kmer_inds[kmer_rc];
        }
        else
        {
            kmer_inds[kmer] = ind;
            kmer_count_len += 1;
            ind += 1;
        }
    }
}

vector<double> count_kmers(string seq)
{
    vector<double> profile(kmer_count_len, 0);
    double total = 0;
    long len = 0;
    u_int64_t val = 0;

    for (int i = 0; i < (int)seq.length(); i++)
    {
        val = (val << 2);
        val = val & ((u_int64_t)pow(4, k_size) - 1);
        val += (seq[i] >> 1 & 3);
        len++;

        if (len == k_size)
        {
            // use val as the kmer for counting
            len--;
            profile[kmer_inds[val]]++;
            total++;
        }
    }

    for (size_t i = 0; i < profile.size(); i++)
    {
        profile[i] /= max(1.0, total);
    }

    return profile;
}

void processLinesBatch(vector<string> &linesBatch, string &outputPath, int threads)
{
    vector<vector<double>> results(linesBatch.size());
    ofstream output;
    output.open(outputPath, ios::out | ios::app);
    string o = "";

#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (size_t i = 0; i < linesBatch.size(); i++)
    {
        results[i] = count_kmers(linesBatch[i]);
    }

    for (size_t i = 0; i < linesBatch.size(); i++)
    {
        for (double d : results[i])
        {
            o += to_string(d);
            o += " ";
        }
        o += "\n";
    }

    output << o;
    results.clear();
    output.close();
}

int main(int argc, char **argv)
{
    vector<string> batch;
    string inputPath, outputPath;
    int threads;
    
    inputPath = argv[1];
    outputPath = argv[2];
    k_size = stoi(argv[3]);
    threads = stoi(argv[4]);

    cout << "INPUT FILE " << inputPath << endl;
    cout << "OUTPUT FILE " << outputPath << endl;
    cout << "K_SIZE " << k_size << endl;
    cout << "THREADS " << threads << endl;

    compute_kmer_inds();

    cout << "Profile Size " << kmer_count_len << endl;
    cout << "Total " << k_size << "-mers " << kmer_inds.size() << endl;

    Seq seq;
    SeqReader reader(inputPath);
    ofstream output;

    output.open(outputPath, ios::out);
    output.close();

    while (reader.get_seq(seq))
    {
        batch.push_back(seq.data);

        if (batch.size() == 10000)
        {
            processLinesBatch(batch, outputPath, threads);
            batch.clear();
        }
    }

    processLinesBatch(batch, outputPath, threads);    
    batch.clear();

    return 0;
}