#!/bin/bash

[ -d mbcclr_utils/bin ] && rm -r mbcclr_utils/bin
mkdir mbcclr_utils/bin

echo "BUILDING THE MetaBCC-LR 15 MER COMPUTATIONS"
g++ mbcclr_utils/search-15mers.cpp -fopenmp -Wall -o mbcclr_utils/bin/search-15mers
echo "BUILDING THE MetaBCC-LR 3 MER COMPUTATIONS"
g++ mbcclr_utils/count-kmers.cpp -Wall -fopenmp -o mbcclr_utils/bin/count-kmers
echo "BUILD FINISHED"