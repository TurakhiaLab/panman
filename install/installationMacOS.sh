#!/bin/bash

# Install dependencies
brew install git cmake wget curl zip unzip tar boost pkg-config protobuf

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
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz

# Run CMake
cmake -DTBB_DIR=${PWD}/oneTBB-2019_U9 \
      -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake \
      -DCMAKE_TOOLCHAIN_FILE=${PWD}/vcpkg/scripts/buildsystems/vcpkg.cmake ..

# Build the project
make -j

# Go back to start directory
cd "$startDir"
