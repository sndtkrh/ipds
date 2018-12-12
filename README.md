# IPDS

## Overview
This is a library for statistical analysis in modern C++.
This library supports following operations:
- basic statistical functions
  - average
  - covariance
  - standard deviation
- operations in linear linear_algebra
  - innner product
  - norm
  - orthonormalization
  - projection
- principal component analysis
- graphically depicting data (svg output support)

## Requirements
To build IPDS, followings are required:
- C++ compiler with C++17 capabilities
  - g++-7 or greater is recomended
- CMake 2.8 or greater
- Boost C++ Libraries
  - uBLAS
Dockerfile is available.

## Build a Docker container
If you don't use Docker, you can skip this step.
```bash
$ docker build --tag ipds .
$ docker run --rm -it --volume ~/:/work --name ipds ipds
```

## Build and Install the Library
```bash
$ git clone git@github.com:sndtkrh/ipds.git
$ cd ./ipds
$ mkdir build
$ cd ./build
$ cmake ..
$ make install
```

## Build and Execute Examples
```bash
$ cd ../examples
$ mkdir build
$ cd ./build
$ cmake ..
$ make
$ ./iris.out < ../iris.txt
```
