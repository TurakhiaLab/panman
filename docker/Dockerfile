FROM --platform=linux/amd64 condaforge/miniforge3:latest
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    apt-add-repository ppa:ubuntu-toolchain-r/test

RUN apt-get update && \
    apt-get install -y gcc-11 \
                    g++-11 git \
                    build-essential bc \
                    pigz \
                    tabix \
                    samtools \
                    make cmake wget curl zip \
                    unzip tar protobuf-compiler \
                    libboost-all-dev pkg-config \
                    git \
                    bash \
                    python3-dev \
                    python3-pip \
                    python3-venv \
                    pybind11-dev \
                    libbz2-dev \
                    libatomic-ops-dev \
                    autoconf \
                    libgsl-dev \
                    zlib1g-dev \
                    libzstd-dev \
                    libjemalloc-dev \
                    libhts-dev \
                    time \
                    pip \
                    libcairo2-dev \
                    parallel \
                    r-base \
                    libclang-dev \
    && apt-get clean \
    && apt-get purge  \
    && rm -rf /var/lib/apt/lists/* && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 100

RUN conda install python=3.10 

WORKDIR /home

# pggb
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 \
    && tar xjf bcftools-1.19.tar.bz2 \
    && cd bcftools-1.19/ && ./configure --prefix=/usr/local/bin/ && make && make install && export PATH=/usr/local/bin/bin:$PATH && cd .. && cp /usr/local/bin/bin/* /usr/local/bin/

RUN git clone --recursive https://github.com/waveygang/wfmash \
   && cd wfmash \
   && git pull \
   && git checkout 517e1bc5c133ecac483a8479c5403f8a13d0fdd5 \
   && git submodule update --init --recursive \
   && sed -i 's/-march=native/-march=sandybridge/g' src/common/wflign/deps/WFA2-lib/Makefile \
   && cmake -H. -DCMAKE_BUILD_TYPE=Generic -DEXTRA_FLAGS='-march=sandybridge -Ofast' -Bbuild && cmake --build build -- -j $(nproc) \
   && cp build/bin/wfmash /usr/local/bin/wfmash \
   && cp build/lib/* /usr/local/lib/ \
   && cd ../

RUN git clone --recursive https://github.com/ekg/seqwish \
   && cd seqwish \
   && git pull \
   && git checkout 0eb6468be0814ab5a0cda10d12aa38cb87d086f1 \
   && git submodule update --init --recursive \
   && cmake -H. -DCMAKE_BUILD_TYPE=Generic -DEXTRA_FLAGS='-march=sandybridge -Ofast' -Bbuild && cmake --build build -- -j $(nproc) \
   && cp bin/seqwish /usr/local/bin/seqwish \
   && cd ../ \
   && rm -rf seqwish

RUN git clone --recursive https://github.com/pangenome/smoothxg \
   && cd smoothxg \
   && git pull \
   && git checkout e93c62356b33c2f6db727452328e09a33bfc82c6 \
   && git submodule update --init --recursive \
   && cmake -H. -DCMAKE_BUILD_TYPE=Generic -DEXTRA_FLAGS='-Ofast' -Bbuild && cmake --build build -- -j $(nproc) \
   && cp bin/smoothxg /usr/local/bin/smoothxg \
   && cp deps/odgi/bin/odgi /usr/local/bin/odgi \
   && cd ../ \
   && rm -rf odgi

# Rust
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN cargo --help

RUN git clone https://github.com/marschall-lab/GFAffix.git \
   && cd GFAffix \
   && git pull \
   && git checkout 7d5fe205a9c46fc4f916174cf8b1ee179e8cdb1f \
   && cargo install --force --path . \
   && mv /root/.cargo/bin/gfaffix /usr/local/bin/gfaffix \
   && cd ../ \
   && rm -rf GFAffix

RUN wget https://github.com/vgteam/vg/releases/download/v1.62.0/vg && chmod +x vg && mv vg /usr/local/bin/vg

RUN git clone https://github.com/pangenome/vcfbub \
   && cd vcfbub \
   && git pull \
   && git checkout db7775f63eab6891acb6000dddfa320146cd7c56 \
   && cargo install --force --path . \
   && mv /root/.cargo/bin/vcfbub /usr/local/bin/vcfbub \
   && cd ../ \
   && rm -rf vcfbub

RUN git clone https://github.com/ekg/fastix.git \
   && cd fastix \
   && git pull \
   && git checkout 331c1159ea16625ee79d1a82522e800c99206834 \
   && cargo install --force --path . && \
   mv /root/.cargo/bin/fastix /usr/local/bin/fastix \
   && cd ../ \
   && rm -rf fastix

RUN wget https://github.com/pangenome/pggb/archive/refs/tags/v0.7.2.zip && unzip v0.7.2.zip && cd pggb-0.7.2 && cp pggb /usr/local/bin && chmod 777 /usr/local/bin

RUN conda create -n panman -y
RUN echo "source activate panman" >> ~/.bashrc
SHELL ["conda", "run", "-n", "panman", "/bin/bash", "-c"]
RUN pip install snakemake

# PanMAN
RUN git clone https://github.com/TurakhiaLab/panman.git && cd panman && mkdir -p build && cd build
RUN curl -O https://capnproto.org/capnproto-c++-1.0.2.tar.gz && tar zxf capnproto-c++-1.0.2.tar.gz && \
   cd capnproto-c++-1.0.2 && sed -i 's/uint64_t traversalLimitInWords = 8 \* 1024 \* 1024;/uint64_t traversalLimitInWords = 8 \* 1024 \* 1024 \* 256;/' $PWD/src/capnp/message.h && ./configure && make -j && make install 

RUN cd panman/build && git clone https://github.com/microsoft/vcpkg.git && apt-get install -y pkg-config && ./vcpkg/bootstrap-vcpkg.sh && ./vcpkg/vcpkg install jsoncpp
RUN cd panman/build && wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz && tar -xvzf 2019_U9.tar.gz
RUN cd panman/build && cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake -DCMAKE_TOOLCHAIN_FILE=${PWD}/vcpkg/scripts/buildsystems/vcpkg.cmake .. && make -j && cp panmanUtils /usr/local/bin && panmanUtils

# PanGraph
RUN wget https://github.com/neherlab/pangraph/archive/refs/tags/0.7.3.zip && unzip 0.7.3.zip && cd pangraph-0.7.3 && rm trace.jl && touch trace.jl && make && export PATH=$PATH:/home/pangraph-0.7.3/pangraph/bin
RUN echo "export PATH=$PATH:/home/pangraph-0.7.3/pangraph/bin" >> ~/.bashrc

RUN conda install bioconda::minimap2
RUN conda install bioconda::mmseqs2
RUN conda install bioconda::mashtree
RUN conda install bioconda::mafft
RUN conda install conda-forge::gdown

# TWILIGHT
RUN apt-get update && \ 
    apt-get install -y wget git build-essential libboost-all-dev libtbb2 cmake protobuf-compiler
RUN git clone https://github.com/TurakhiaLab/TWILIGHT.git -b experiment/sars && \
    mkdir ./TWILIGHT/build && \
    cd ./TWILIGHT/build && \ 
    wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz && \
    tar -xvzf 2019_U9.tar.gz && \
    cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  .. && \
    make -j && \
    cp twilight /usr/local/bin

RUN echo "export LD_LIBRARY_PATH=/home/TWILIGHT/build/tbb_cmake_build/tbb_cmake_build_subdir_release/:$LD_LIBRARY_PATH" >> ~/.bashrc

RUN wget https://mol.ax/rxiv/3seq_build_170612.zip && \
    unzip 3seq_build_170612.zip && \
    mv '3seq build 170612' 3Seq && \
    cd 3Seq && make -j && cp 3seq /usr/local/bin && cd ../ && rm -r 3*
   