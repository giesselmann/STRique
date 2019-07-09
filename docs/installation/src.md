# Source

STRique can be installed without root privileges as it's maintaining most of its dependencies by building required libraries from source in a user accessible path. Currently the following list of system wide packages is required to build all included software packages (names on MacOS might be different):

* git wget
* gcc g++ make cmake
* python3-dev

These packages are likely present in most production environments. Please also refer to the Dockerfile of the package repository.
