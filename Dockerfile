# GCC 7
FROM ubuntu:16.04

RUN apt-get update \
  && apt-get install --yes --no-install-recommends \
    build-essential \
    cmake \
    emacs \
    git \
    libboost-all-dev \
    software-properties-common \
    ssh \
  && rm -rf /var/lib/apt/lists/*

RUN add-apt-repository ppa:ubuntu-toolchain-r/test

RUN apt-get update \
  && apt-get install --yes --no-install-recommends \
    gcc-7 \
    g++-7 \
  && rm -rf /var/lib/apt/lists/*

ENV CC=/usr/bin/gcc-7
ENV CXX=/usr/bin/g++-7
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 70 --slave /usr/bin/g++ g++ /usr/bin/g++-7

RUN git clone http://git.tiker.net/trees/boost-numeric-bindings.git \
  && cd boost-numeric-bindings \
  && ./configure --prefix=/usr/local/ \
  && make install \
  && cd ../ \
  && rm -rf boost-numeric-bindings

ENV CPLUS_INCLUDE_PATH="$CPLUS_INCLUDE_PATH:/usr/local/include/boost-numeric-bindings"

RUN apt-get update \
  && apt-get install --yes --no-install-recommends \
    liblapack3 \
    liblapack-dev \
    libopenblas-base \
    libopenblas-dev \
    liblapacke-dev \
    liblapack-dev \
  && rm -rf /var/lib/apt/lists/*
