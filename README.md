# General Information

http://www.fesslix.org (offical website)


# Copying

Fesslix is distributed under the terms of the General Public License which
is compliant with the guidelines of the Open Source and Free Software
Fundations. See the file COPYING for details.


# Installation

http://www.fesslix.org/installation


# Building the package

## Environment & Dependencies

Install the package 'docker'.

python -m venv .venv-fesslix
source .venv-fesslix/bin/activate

pip install build scikit-build-core pybind11
pip install twine cibuildwheel

## Local build

python -m build
pip install .

## Precompiled Binaries
https://pypi.org/account/register/

python -m cibuildwheel --output-dir dist
twine upload dist/*

