#!/bin/bash
gcc --version
g++ --version
[ ! -d vcpkg ] && git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.sh
./vcpkg install gsl
find ./ -name libgsl*
pwd
