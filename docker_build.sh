#!/bin/bash

echo "Building only the C++ components"

[ -d mbcclr_utils/bin ] && rm -r mbcclr_utils/bin
mkdir mbcclr_utils/bin  

echo "Linux Build"
# Linux
echo "BUILDING THE MetaBCC-LR 15 MER COMPUTATIONS"
g++ mbcclr_utils/search-15mers.cpp -fopenmp -lpthread -Wall -o mbcclr_utils/bin/search-15mers -lz -O3
g++ mbcclr_utils/count-15mers.cpp -fopenmp -lpthread -Wall -o mbcclr_utils/bin/count-15mers -lz -O3
echo "BUILDING THE MetaBCC-LR 3 MER COMPUTATIONS"
g++ mbcclr_utils/count-kmers.cpp -Wall -fopenmp -lpthread -o mbcclr_utils/bin/count-kmers -lz -O3
echo "BUILD FINISHED"
