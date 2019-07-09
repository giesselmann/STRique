# Docker

The all-in-one Docker image of STRique contains all its dependencies and is build automatically on [dockerHub](https://hub.docker.com/r/giesselmann/STRique). The compressed docker image size is ~700 MB. From the Docker shell run:

    docker pull giesselmann/strique

The container is tagged with the package version. To obtain a specific version run e.g.:

    docker pull giesselmann/strique:v0.3.0

Docker images can get access to the hosts file system with mounts of type *bind*. Test the container by mounting your current directory:

    docker run -it --mount type=bind,source=$(pwd),target=/host giesselmann/strique:v0.3.0

Inside of the container type

    ls -l /host

to see the files of the host system from where the container was started. Leave a running container with *exit*.
