FROM --platform=linux/amd64 condaforge/miniforge3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && \
    apt-get install -y gcc-11 g++-11 git build-essential \
                   cmake wget curl zip \
                   unzip tar protobuf-compiler \
                   libboost-all-dev pkg-config && \
    apt-get clean
    # update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100 && \
    # update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 100


WORKDIR /home

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

# Minimap2
RUN conda install bioconda::minimap2

# MMSEQ
RUN conda install bioconda::mmseqs2

# MashTree
RUN conda install bioconda::mashtree

# WFMash
RUN conda install bioconda::wfmash

# WFMash
RUN conda install bioconda::mafft

# PGGB and VG
RUN conda install bioconda::pggb

RUN conda install conda-forge::gdown
