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


# Commiting code

git push -u origin main


