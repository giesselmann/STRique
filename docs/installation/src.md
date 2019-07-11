# Source

STRique can be installed without root privileges as it's maintaining its dependencies by building required libraries from source in a user accessible path. Currently the following list of system wide packages is required to build all included software packages (names on macOS might be different):

* git wget
* gcc g++ make cmake
* python3-dev

These packages are likely present in most production environments. Please also refer to the Dockerfile of the package repository for minimal system requirements.

??? info "Python3 macOS"
    While git, gcc and make are present on macOS, cmake and python3 might not be. If not already done consider  installing the [brew](https://wsvincent.com/install-python3-mac/) package management and run ```brew install  cmake python3``` afterwards. It can be necessary to modify your PATH variable to find the python3 interpreter.  Before continuing with this guide make sure that either ```python --version``` or ```python3 --version```   yields the expected result of 3.x.x

## Dependencies
**Python 3.5 or higher**

    pomegranate
    numpy, scipy, scikit-image
    h5py

**C++**

    g++ >= 6 (C++14 support)
    cmake
    SeqAn2
    Pybind11

Dependencies get downloaded and build by the setup script.

## Installation

Please consider to install STRique into a separate virtual environment to avoid conflicting dependencies with other packages.

Create a virtual python environment:

```
python3 -m venv /path/to/your/environment
cd /path/to/your/environment
source bin/activate
```

In order to download, build and install STRique, execute the following commands:

```
mkdir -p ~/src && cd ~/src
git clone --recursive https://github.com/giesselmann/STRique
cd STRique
pip3 install -r requirements.txt
python3 setup.py install
```

If the process didn't report any errors, continue in the [test](test.md) section to verify the installation.
