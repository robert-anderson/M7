# M7
Many-body Stochastic Expectation Value Estimation Networks (M7) is a stochastic electronic structure program written in C++ which primarily implements the Full Configuration Interaction Quantum Monte Carlo (FCIQMC) method.

# Installation
To install M7, you will need:
* C++11 compiler which implements MPI
* CMake installation at least as recent as version 3.10
* LAPACK implementation

Retrieve the program code from the VCS:
```bash
git clone https://github.com/robert-anderson/M7.git
```
Then create a directory in which CMake can build its targets---typically under the repository root directory---and point CMake to the top level CMakeLists.txt file of the project.
```bash
cd M7
mkdir build
cd build
cmake ..
make -j TARGET
```
Where TARGET is either debug, release or unittest.
The binary will then appear at either src/debug, src/release, or test/unittest.

Parallel HDF5 is a non-prerequisite dependency which is downloaded, compiled, and linked as a static library during the build process.
If using the Intel C++ compiler, make sure the compilevars.sh has been sourced to correctly configure the library path environment variables prior to invoking CMake.
Otherwise, HDF5 may fail to find Intel symbols such as `_intel_fast_memcpy` at link time.

Also, the HDF5 build script may quit reporting that MPI_C_FOUND is not set, meaning that CMake's FIND_PACKAGE feature could not determine the path of the MPI library.
Ensuring that the MPI library root is included in the PATH environment variable should be sufficient to overcome this.
This issue may arise in IDEs where the MPI compiler executable has been configured as the default compiler, which allows the M7 CMakeLists.txt file to find the MPI path, but fails for the HDF5 external project.
In such cases, set the PATH variable within the IDE's toolchain dialogue to contain the MPI library root.


# Usage
The program takes all run time variables through the command line interface.
The command line options available can be listed by executing the debug or release binaries with the -h or --help option.

Complex configurations can become cumbersome to manage on the command line directly, so particular run time configurations can easily be split over multiple lines in a shell script.

If the user chooses, the python/ directory contains a simple and convenient wrapper which encapsulates all options, and can handle the execution of distributed instances of M7.
