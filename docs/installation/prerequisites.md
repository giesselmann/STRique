# Installation

*Important*: If you encounter any trouble within the installation, open an issue on [github](https://github.com/giesselmann/STRique/issues) to make others aware of the problem or write an E-Mail to ```giesselmann[at]molgen.mpg.de```. There is nothing more annoying than software that does not work.

In general STRique is python based and therefore available on Windows, MacOS and Linux. However the included signal alignment is a C++ extension which needs to be build on the target system. We support and test the fully featured package on Linux and MacOS and through containerization on Windows.

<center>

| Method   | Source 	| Docker 	|
|---------|:------:	|:------:	|
| Linux   |   yes  	|   yes  	|
| MacOS   |   yes  	|   yes  	|
| Windows |   no   	|   yes  	|

</center>

**Recommendation:** If possible use the Docker version of STRique. Singularity might be an alternative in case the root access required by Docker is not possible on your machine. Lastly with the source build you can integrate STRique in any (virtual) python3 environment.
