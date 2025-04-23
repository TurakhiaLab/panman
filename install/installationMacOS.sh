#!/bin/bash

# Install dependencies
brew install git cmake wget curl zip unzip boost pkg-config protobuf libtool automake autoconf nasm tbb

# wget https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz
# tar -xvf v2.30.0.tar.gz
# pushd isa-l-2.30.0
# ./autogen.sh
# ./configure --prefix=$(brew --prefix) --libdir=$(brew --prefix)/lib
# make -j2
# make install
# popd

# Set start directory
startDir=$(pwd)
cd "$(dirname "$0")"
mkdir -p ../build
cd ../build

# Install capnp manually if newer version needed
curl -O https://capnproto.org/capnproto-c++-1.0.2.tar.gz
tar zxf capnproto-c++-1.0.2.tar.gz
cd capnproto-c++-1.0.2
sed -i '' 's/uint64_t traversalLimitInWords = 8 \* 1024 \* 1024;/uint64_t traversalLimitInWords = 8 \* 1024 \* 1024 \* 256;/' src/capnp/message.h
./configure
make -j
sudo make install
which capnp
cd ../

# Clone and bootstrap vcpkg
git clone https://github.com/microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install jsoncpp
./vcpkg/vcpkg install abseil

# Download and extract oneTBB
# wget https://github.com/oneapi-src/oneTBB/releases/download/2019_U9/tbb2019_20191006oss_mac.tgz
# tar -xvzf tbb2019_20191006oss_mac.tgz
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xzf 2019_U9.tar.gz
cd oneTBB-2019_U9
make -j

sudo cp -r include/tbb /usr/local/include/
sudo cp build/*_release/*.dylib /usr/local/lib/

# Run CMake
# cmake -DTBB_DIR=${PWD}/tbb2019_20191006oss \
#       -DCMAKE_PREFIX_PATH=${PWD}/tbb2019_20191006oss/cmake \
cmake -DCMAKE_TOOLCHAIN_FILE=${PWD}/vcpkg/scripts/buildsystems/vcpkg.cmake ..

# Build the project
make -j

# Go back to start directory
cd "$startDir"
