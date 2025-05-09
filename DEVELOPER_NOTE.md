# Building the package

## Environment & Dependencies

TODO: needed? Build is documented on: fesslix-docu/docs/_build/html/start_installation.html#install-as-a-local-build

Install the package 'docker'.

python -m venv .venv-fesslix
source .venv-fesslix/bin/activate

pip install build scikit-build-core pybind11
pip install twine cibuildwheel

## fast install with KDevelop (for development)

set CMAKE_INSTALL_PREFIX to e.g., ~/.venv-fesslix/lib/python3.13/site-packages
then install with SHIFT+F8


## Precompiled Binaries
https://pypi.org/account/register/

python -m cibuildwheel --output-dir dist
twine upload dist/*

## GitHub Actions

### cibuildwheel 
https://cibuildwheel.pypa.io/en/stable/faq/
https://cibuildwheel.pypa.io/en/stable/working-examples/

### Trigger the workflow on push and pull request events to the main branch

```
on:
  push:
    branches:
      - main
  pull_request:
    branches:
      - main
```

### Trigger the workflow action whenever a new version tag is committed

```
on:
  push:
    tags:
      - 'v*'  # Trigger only on version tags
```

To create a tag:
```
git tag v1.0.0  # Create a new version tag

```
To transfer the tag to the server:
```
git push origin v1.5
```
To transfer all tags to the server:
```
git push origin --tags
```



# CMake
Output messages in CMake:

```
message("202502240926 Boost_INCLUDE_DIR=${Boost_INCLUDE_DIR}")
```



# Commiting code

git push -u origin main

## commit messages

increase patch number
attempt package build


# Check for memory leaks




# Documentation

see <https://github.com/fesslix/fesslix-docu>

