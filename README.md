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

## Build and Install the Library
```bash
$ git clone git@github.com:sndtkrh/ipds.git
$ cd ./ipds
$ mkdir build
$ cd ./build
$ cmake ..
$ make
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
