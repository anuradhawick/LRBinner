#!/bin/bash

[ -d mbcclr_utils/bin ] && rm -r mbcclr_utils/bin
mkdir mbcclr_utils/bin    

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
    g++ mbcclr_utils/search-15mers.cpp -fopenmp -lpthread -Wall -o mbcclr_utils/bin/search-15mers -lz -O3
    g++ mbcclr_utils/count-15mers.cpp -fopenmp -lpthread -Wall -o mbcclr_utils/bin/count-15mers -lz -O3
    echo "BUILDING THE MetaBCC-LR 3 MER COMPUTATIONS"
    g++ mbcclr_utils/count-kmers.cpp -Wall -fopenmp -lpthread -o mbcclr_utils/bin/count-kmers -lz -O3
    echo "BUILD FINISHED"
    ;;
esac