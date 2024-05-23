# Install dependencies

# Build
startDir=$pwd
cd $(dirname "$0")
mkdir -p ../build
cd ../build
git clone https://github.com/microsoft/vcpkg.git
./vcpkg/vcpkg install jsoncpp
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake -DProtobuf_PROTOC_EXECUTABLE=/usr/local/bin/protoc -DCMAKE_TOOLCHAIN_FILE=${PWD}/vcpkg/scripts/buildsystems/vcpkg.cmake ..
make -j

cd $startDir