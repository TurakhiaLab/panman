# brew install gcc cmake boost abseil

rm -rf build && mkdir build

cd build

    wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
    wget https://github.com/protocolbuffers/protobuf/releases/download/v25.1/protobuf-25.1.tar.gz

    tar xvzf 2019_U9.tar.gz
    tar xvzf protobuf-25.1.tar.gz

    cd ${PWD}/protobuf-25.1
        rm -rf build && mkdir build
        cd build
            cmake -DBUILD_SHARED_LIBS=ON -Dprotobuf_BUILD_TESTS=OFF -Dprotobuf_ABSL_PROVIDER=package ..
            make -j8
        cd ..                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
    cd ..

    cmake -DCMAKE_PREFIX_PATH="${PWD}/oneTBB-2019_U9/cmake" \
        -DProtobuf_DIR="${PWD}/protobuf-25.1/cmake/" \
        -DTBB_DIR=${PWD}/oneTBB-2019_U9 ..
    
    make -j8