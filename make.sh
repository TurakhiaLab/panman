#brew install gcc cmake boost abseil

rm -rf build && mkdir build

cd build

    wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
    tar xvzf 2019_U9.tar.gz
    git clone https://github.com/rvaser/spoa && cd spoa
    cmake -B build -DCMAKE_BUILD_TYPE=Release -Dspoa_install=ON
    make -C build && cd ..

    cmake -DTBB_DIR=${PWD}/oneTBB-2019_U9 -Dspoa_DIR=${PWD}/spoa/  ..
    
    make -j8