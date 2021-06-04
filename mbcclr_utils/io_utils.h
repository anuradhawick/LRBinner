#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

class Seq
{
public:
    string id;
    string data;
};

// class SeqReader
// {
// private:
//     string path;
//     ifstream file_handler;
//     string seq_id = "";
//     string seq_data = "";
//     string line;
//     int line_no = -1;
//     int mode;

// public:
//     SeqReader(string path)
//     {
//         this->path = path;
//         file_handler.open(path, ios::in);
//     }

//     ~SeqReader()
//     {
//         file_handler.close();
//     }

//     bool get_seq(Seq &seq)
//     {

//         while (getline(file_handler, line))
//         {
//             if (line.length()==0)
//             {
//                 continue;
//             }
//             line_no++;
//             // detect format
//             if (line_no == 0)
//             {
//                 if (line[0] == '>')
//                 {
//                     mode = 1;
//                 }
//                 else
//                 {
//                     mode = 2;
//                 }
//             }
//             // fasta
//             if (mode == 1)
//             {
//                 // seeing a new entry
//                 if (line[0] == '>')
//                 {
//                     if (seq_id.length() > 0)
//                     {
//                         // prepare to return
//                         seq.id = seq_id;
//                         seq.data = seq_data;

//                         // store the observed
//                         seq_id = line.substr(1);
//                         seq_data = "";

//                         return true;
//                     }
//                     else
//                     {
//                         seq_id = line.substr(1);
//                     }
//                 }
//                 else
//                 {
//                     seq_data += line;
//                 }
//             }
//             // fastq
//             else if (mode == 2)
//             {
//                 // new entry seen
//                 if (line_no % 4 == 0)
//                 {
//                     seq_id = line.substr(1);
//                     seq_data = "";
//                 }
//                 // new entry end
//                 else if (line_no % 4 == 1)
//                 {
//                     // prepare to return
//                     seq.id = seq_id;
//                     seq.data = line;

//                     return true;
//                 }
//             }
//         }

//         if (mode == 2)
//         {
//             return false;
//         }
//         else if (seq_id.length() > 0)
//         {
//             seq.id = seq_id;
//             seq.data = seq_data;

//             seq_id = "";

//             return true;
//         }
//         else
//         {
//             return false;
//         }
//     }
// };

class SeqReader
{
private:
    gzFile fp;
    kseq_t *ks;
    int ret;

public:
    SeqReader(string path)
    {
        fp = gzopen(path.c_str(), "r");
        ks = kseq_init(fp);
    }

    ~SeqReader()
    {
        kseq_destroy(ks);
        gzclose(fp);
    }

    bool get_seq(Seq &seq)
    {
        if ((ret = kseq_read(ks)) >= 0)
        {
            // using raw pointers can have unwated side effects
            // better to copy than debug
            seq.data = string(ks->seq.s);
            seq.id = string(ks->name.s);
            return true;
        }
        return false;
    }
};