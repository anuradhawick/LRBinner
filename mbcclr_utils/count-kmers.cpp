#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <queue>
#include <mutex>
#include <thread>
#include <condition_variable>
#include "io_utils.h"

using namespace std;

uint32_t kmer_count_len = 0;
map<uint32_t, uint32_t> kmer_inds;
uint32_t *kmer_inds_index;
uint32_t k_size;
queue<string> reads_queue;
mutex mux;
condition_variable condition;
volatile bool terminate_threads;

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

    kmer_inds_index = new uint32_t[kmer_inds.size()];

    for (auto its = kmer_inds.begin(); its != kmer_inds.end(); its++)
    {
        kmer_inds_index[its->first] = its->second;
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
            profile[kmer_inds_index[val]]++;
            total++;
        }
    }

    for (size_t i = 0; i < profile.size(); i++)
    {
        profile[i] /= max(1.0, total);
    }

    return profile;
}

void processLinesBatch(vector<string> &linesBatch, string &output_path, int threads)
{
    vector<vector<double>> results(linesBatch.size());
    ofstream output;
    output.open(output_path, ios::out | ios::app);
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

void off_load_process(string &output, int &threads)
{
    string seq;
    vector<string> batch;

    while (true)
    {
        {
            unique_lock<mutex> lock(mux);

            while (reads_queue.size() > 0)
            {
                seq = reads_queue.front();
                batch.push_back(seq);
                reads_queue.pop();

                if (batch.size() == 10000)
                {
                    break;
                }
            }
        }

        condition.notify_all();

        if (batch.size() > 0)
        {
            processLinesBatch(batch, output, threads);
            batch.clear();
        }

        {
            unique_lock<mutex> lock(mux);
            if (terminate_threads && reads_queue.size() == 0)
            {
                break;
            }
        }
    }
}

void io_thread(string &file_path)
{
    SeqReader reader(file_path);
    Seq seq;
    int count = 0;

    while (reader.get_seq(seq))
    {
        {
            unique_lock<mutex> lock(mux);
            condition.wait(lock, [] { return reads_queue.size() < 50000; });
            reads_queue.push(seq.data);
        }
        count++;

        cout << "Loaded Reads " << count << "       \r" << flush;
    }

    cout << endl;

    terminate_threads = true;
}

int main(int argc, char **argv)
{
    vector<string> batch;
    string input_path, output_path;
    int threads;

    input_path = argv[1];
    output_path = argv[2];
    k_size = stoi(argv[3]);
    threads = stoi(argv[4]);

    cout << "INPUT FILE " << input_path << endl;
    cout << "OUTPUT FILE " << output_path << endl;
    cout << "K_SIZE " << k_size << endl;
    cout << "THREADS " << threads << endl;

    compute_kmer_inds();

    cout << "Profile Size " << kmer_count_len << endl;
    cout << "Total " << k_size << "-mers " << kmer_inds.size() << endl;

    ofstream output(output_path, ios::out);

    thread iot(io_thread, ref(input_path));
    thread process(off_load_process, ref(output_path), ref(threads));

    iot.join();
    process.join();

    return 0;
}