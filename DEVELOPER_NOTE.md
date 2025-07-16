# Building the package

## Environment & Dependencies

TODO: needed? Build is documented on: fesslix-docu/docs/_build/html/start_installation.html#install-as-a-local-build

Install the package 'docker'.

python -m venv .venv-fesslix
source .venv-fesslix/bin/activate

OBSOLETE?!
#pip install build scikit-build-core pybind11
#pip install twine cibuildwheel

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
To list all tags:
```
git tag
```
To list tags that match a pattern:
```
git tag -l "v*"
```


# CMake
Output messages in CMake:

```
message("202502240926 Boost_INCLUDE_DIR=${Boost_INCLUDE_DIR}")
```



# Commiting code

git push -u origin main

## common commit messages

versioning: increase patch number
versioning: increase minor number
versioning: increase major number
attempt package build
change CI-build trigger from 'push' to 'tag'

## work on features in a separate branch

### list existing branches

locally:
    git branch
    
remotely
    git branch -r
    

### Create a new branch

git switch -c GPR

### Switch branches

git switch main
git switch GPR

### Save local changes without commiting

git stash

»»» now, you can switch branches and do other stuff (on the version of the latest commit)
»»» afterwards, make sure you are in the intendend branch

git stash pop


### merge changes in the main branch into the feature branch

git switch GPR

#### merge from remote main branch
git fetch origin
git merge origin/main

#### merge from local main branch
git merge main


### merge changes back to main branch

git switch main
git pull origin main
git merge GPR
git branch -d GPR



# Check for memory leaks




# Documentation

see <https://github.com/fesslix/fesslix-docu>

