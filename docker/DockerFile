FROM ubuntu:20.04

RUN apt update 
RUN apt install -y git build-essential cmake wget curl zip unzip tar protobuf-compiler libboost-all-dev pkg-config capnproto

WORKDIR /HOME

RUN git clone https://github.com/TurakhiaLab/panman.git
RUN panman/install/installationUbuntu.sh
