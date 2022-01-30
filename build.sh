#!/bin/bash

cpath=$(pwd)
[ -d mbcclr_utils/bin ] && rm -r mbcclr_utils/bin
mkdir mbcclr_utils/bin    

[ -d auxiliary/FragGeneScan1.31 ] && rm -r auxiliary/FragGeneScan1.31

wget -O auxiliary/FragGeneScan1.31.tar.gz https://sourceforge.net/projects/fraggenescan/files/FragGeneScan1.31.tar.gz
tar -xzf auxiliary/FragGeneScan1.31.tar.gz -C auxiliary/ && rm auxiliary/FragGeneScan1.31.tar.gz
cd auxiliary/FragGeneScan1.31 && make clean && make fgs && cd ../..
cd $cpath

[ -d auxiliary/hmmer-3.3.2 ] && rm -r auxiliary/hmmer-3.3.2

wget -O auxiliary/hmmer.tar.gz http://eddylab.org/software/hmmer/hmmer.tar.gz
tar -xzf auxiliary/hmmer.tar.gz -C auxiliary/ && rm auxiliary/hmmer.tar.gz
cd auxiliary/hmmer-3.3.2 && ./configure && make && cd ../..
cd $cpath

case $1 in

  osx | macos)
    echo "OSX Build"
    # OSX Build (modify the include path to suit you; you might need to run brew install libomp or brew install llvm)
    echo "BUILDING THE MetaBCC-LR 15 MER COMPUTATIONS"
    clang++ mbcclr_utils/search-15mers.cpp -lomp -fopenmp -lpthread -Wall -o mbcclr_utils/bin/search-15mers -I/usr/local/include -L/usr/local/lib -lz -O3
    clang++ mbcclr_utils/count-15mers.cpp -lomp -fopenmp -lpthread -Wall -o mbcclr_utils/bin/count-15mers -I/usr/local/include -L/usr/local/lib -lz -O3
    echo "BUILDING THE MetaBCC-LR 3 MER COMPUTATIONS"
    clang++ mbcclr_utils/count-kmers.cpp -Wall -lomp -fopenmp -lpthread -o mbcclr_utils/bin/count-kmers -I/usr/local/include -L/usr/local/lib -lz -O3
    echo "BUILD FINISHED"
    ;;

  *)
    echo "Linux Build"
    # Linux
    echo "BUILDING THE MetaBCC-LR 15 MER COMPUTATIONS"
    g++ mbcclr_utils/search-15mers.cpp -fopenmp -lpthread -Wall -o mbcclr_utils/bin/search-15mers -lz -O3 -std=c++11
    g++ mbcclr_utils/count-15mers.cpp -fopenmp -lpthread -Wall -o mbcclr_utils/bin/count-15mers -lz -O3 -std=c++11
    echo "BUILDING THE MetaBCC-LR 3 MER COMPUTATIONS"
    g++ mbcclr_utils/count-kmers.cpp -Wall -fopenmp -lpthread -o mbcclr_utils/bin/count-kmers -lz -O3 -std=c++11
    echo "BUILD FINISHED"
    ;;
esac
