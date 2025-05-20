#!/bin/bash
ls /opt/python/
python --version
gcc --version
g++ --version
[ ! -d vcpkg ] && git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.sh
./vcpkg install gsl
./vcpkg install boost-format boost-math boost-concept-check boost-random boost-algorithm
#./vcpkg install nlopt
find ./ -name libgsl*
pwd
cd ..
