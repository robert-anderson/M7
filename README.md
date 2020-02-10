# M7
Many-body Stochastic Expectation Value Estimation Networks (M7) is a stochastic electronic structure program written in C++ which primarily implements the Full Configuration Interaction Quantum Monte Carlo (FCIQMC) method.

# Installation
To install M7, you only need:
* git
* A C++11 compiler which implements MPI
* A CMake installation at least as recent as version 3.14

The program depends on several VCS-hosted packages, which must be collected through git with the following clone command:
```bash
git clone --recurse-submodules https://github.com/robert-anderson/M7.git
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

# Usage
The program takes all run time variables through the command line interface.
The command line options available can be listed by executing the debug or release binaries with the -h or --help option.

Complex configurations can become cumbersome to manage on the command line directly, so particular run time configurations can easily be split over multiple lines in a shell script.

If the user chooses, the python/ directory contains a simple and convenient wrapper which encapsulates all options, and can handle the execution of distributed instances of M7.
