#!/bin/bash

# Install dependencies
brew install git cmake wget curl zip unzip tar boost pkg-config protobuf rsync openmpi libtool automake autoconf nasm

wget https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz
tar -xvf v2.30.0.tar.gz
pushd isa-l-2.30.0
./autogen.sh
./configure --prefix=$(brew --prefix) --libdir=$(brew --prefix)/lib
make -j2
make install
popd

# Set start directory
startDir=$(pwd)
cd "$(dirname "$0")"
mkdir -p ../build
cd ../build

# Install capnp manually if newer version needed
curl -O https://capnproto.org/capnproto-c++-1.0.2.tar.gz
tar zxf capnproto-c++-1.0.2.tar.gz
cd capnproto-c++-1.0.2
./configure
make -j6 check
sudo make install
which capnp
cd ../

# Clone and bootstrap vcpkg
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.sh
cd ..

# Install jsoncpp via vcpkg
./vcpkg/vcpkg install jsoncpp

# Download and extract oneTBB
wget https://github.com/oneapi-src/oneTBB/releases/download/2019_U9/tbb2019_20191006oss_mac.tgz
tar -xvzf tbb2019_20191006oss_mac.tgz

# Run CMake
cmake -DTBB_DIR=${PWD}/tbb2019_20191006oss \
      -DCMAKE_PREFIX_PATH=${PWD}/tbb2019_20191006oss/cmake \
      -DCMAKE_TOOLCHAIN_FILE=${PWD}/vcpkg/scripts/buildsystems/vcpkg.cmake ..

# Build the project
make -j

# Go back to start directory
cd "$startDir"
