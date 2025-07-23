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
cd ../..
echo $PWD
echo "========================="
echo "TBB » build from source"
echo "========================="
git clone https://github.com/oneapi-src/oneTBB.git
cd oneTBB
mkdir build && cd build
cmake ..  -DTBB_TEST=OFF -DTBB_STRICT=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=$PWD/../install -DCMAKE_CXX_FLAGS="-Wno-error"
make -j
make install
cd ../..
echo $PWD
