#brew install gcc cmake boost abseil


cd build

    cmake -DTBB_DIR=${PWD}/oneTBB-2019_U9 -Dspoa_DIR=${PWD}/spoa/  ..
    
    make -j8
    