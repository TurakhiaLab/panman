conda activate panmat
mamba install -c conda-forge abseil-cpp
mamba install -c anaconda boost

rm -rf build && mkdir build && cd build

wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
wget https://github.com/protocolbuffers/protobuf/releases/download/v25.1/protobuf-25.1.tar.gz

tar xvzf 2019_U9.tar.gz
tar xvzf protobuf-25.1.tar.gz

cmake -DCMAKE_PREFIX_PATH="${PWD}/oneTBB-2019_U9/cmake" \
    -DProtobuf_DIR="${PWD}/protobuf-25.1/cmake/" \
    -DTBB_DIR=${PWD}/oneTBB-2019_U9 ..