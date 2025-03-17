# Install dependencies

sudo apt install -y git build-essential cmake wget curl zip unzip tar libboost-all-dev pkg-config protobuf-compiler


# Build
startDir=$pwd
cd $(dirname "$0")
mkdir -p ../build
cd ../build

# install capnp
curl -O https://capnproto.org/capnproto-c++-1.0.2.tar.gz
tar zxf capnproto-c++-1.0.2.tar.gz
cd capnproto-c++-1.0.2
./configure
make -j6 check
sudo make install
which capnp
cd ../

git clone https://github.com/microsoft/vcpkg.git
sudo apt-get install pkg-config
./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install jsoncpp

wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz

cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake -DCMAKE_TOOLCHAIN_FILE=${PWD}/vcpkg/scripts/buildsystems/vcpkg.cmake ..

make -j

cd $startDir
