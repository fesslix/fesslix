#!/bin/bash
echo "========================="
echo "NLopt » build from source"
echo "========================="
git clone https://github.com/stevengj/nlopt.git
cd nlopt
mkdir build && cd build
cmake .. -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=$PWD/../install
make -j
make install
echo $PWD
