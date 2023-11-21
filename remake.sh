
cd build

cmake -DCMAKE_PREFIX_PATH="${PWD}/oneTBB-2019_U9/cmake" \
    -DProtobuf_DIR="${PWD}/protobuf-25.1/cmake/" \
    -DTBB_DIR=${PWD}/oneTBB-2019_U9 ..

make -j8