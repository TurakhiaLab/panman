# Installation Methods

## Using installation script (requires sudo access)

0. Dependencies
    i. Git

1. Clone the repository
```bash
git clone https://github.com/TurakhiaLab/panman.git
cd panman
```
2. Run the installation script
```bash
chmod +x install/installationUbuntu.sh
./install/installationUbuntu.sh
```
3. Run panmanUtils
```bash
cd build
./panmanUtils --help
```
!!!Note
    <i>panmanUtils</i> is built using CMake and depends upon libraries such as Boost, cap'n proto, etc, which are also installed in `installationUbuntu.sh`. If users face version issues, try using the docker methods detailed below.

## Using Docker Image

To use <i>panmanUtils</i> in a docker container, users can create a docker container from a docker image, by following these steps

0. Dependencies
    i. Docker
1. Pull the PanMAN docker image from DockerHub
```bash
docker pull swalia14/panman:latest
```
2. Build and run the docker container
```bash
docker run -it swalia14/panman:latest
```
3. Run panmanUtils
```bash
# Insider docker container
cd /home/panman/build
./panmanUtils --help
```
!!!Note
 The docker image comes with preinstalled panmanUtils and other tools such as PanGraph, PGGB, and RIVET.

## Using DockerFile
Docker container with preinstalled <i>panmanUtils</i> can also be built from DockerFile by following these steps

0. Dependencies
    i. Docker
    ii. Git
1. Clone the repository
```bash
git clone https://github.com/TurakhiaLab/panman.git
cd panman
```
2. Build a docker image
```bash
cd docker
docker build -t panman .
```
3. Build and run docker container
```bash
docker run -it panman
```
4. Run panmanUtils
```bash
# Insider docker container
cd /home/panman/build
./panmanUtils --help
```
