# M7
Many-body Stochastic Expectation Value Estimation Networks (M7) is a stochastic electronic structure program written in C++ which primarily implements the Full Configuration Interaction Quantum Monte Carlo (FCIQMC) method.

# Features 
* MPI parallelization with dynamic load balancing
* YAML-based configuration syntax
* Complex-valued wavefunctions and Hamiltonians can be enabled at compile time
* Compatibility with standard FCIDUMPs in plain text and MOLCAS HDF5 formats
* Sz-unrestricted integrals can be specified in either spin-major or spin-minor indexing schemes, or the Molpro spin blocks protocol
* Sz non-conserving Hamiltonians supported e.g. SOC or Dirac-Coulomb
* All combinations of 2-electron integral permutational symmetries (real orbital, hermiticity, dummy variable) are supported and automatically detected from the FCIDUMP
* Pre-computed heatbath excitation generators
* Initiator adaptation
* In addition to normal fermionic FCIQMC, users can compile with bosons enabled to perform FCIQMC in a space of determinant-permanent products 
* HDF5 interface for saving/loading walker distributions and RDMs
* spin-resolved and spinfree RDMs for fermion ranks 1, 2, 3, and MRPT2 intermediates, as well as multidimensional expectation values for fermion-boson product operators
* Refreshable semi-stochastic spaces

# Building the binaries
M7's build system employs CMake. The compile time variables are:
* `CMAKE_BUILD_TYPE` (default `Release`) either `Release` or `Debug`. `Debug` turns off all compiler optimizations and enabled execution of all macros beginning with `DEBUG_ASSERT_`. `Release` does not include execution of these checks and enables the highest level of compiler optimizations. 
* `BUILD_HDF5` (default `off`) if this is set to `on`, HDF5 will be downloaded and built as part of the M7 build. Otherwise, CMake's `find_package` will look for a suitable extant installation. This installation must be compiled with parallel I/O support or `M7_lib` will not compile
* `MBF_TYPE` (default `fermion`) MBF is used throughout the program's source code and outputs to refer to a Many-body Basis Function, a general term which includes Slater determinants, bosonic permanants, and determinant-permanent products. These are selected at compile time with `fermion`, `boson`, and `fermion-boson` respectively
* `ENABLE_COMPLEX` (default `off`) if this is set to `on`, the walker populations and integral storage vectors will gain an imaginary part. This is required if the user's FCIDUMP is complex-valued.
* `ENABLE_LOCAL_LOGGING` (default `on` if `CMAKE_BUILD_TYPE`=`Debug`, else `off`) if this is set to `on`, each MPI rank produces its own local logs to `M7.*.log` files, which are primarily useful in debugging contexts. In any case, the rank-reduced logs emitted by the root MPI rank are always output to stdout with ANSI formatting, and reproduced without formatting in `M7.log`

The build first produces a statically linked library `M7_lib.a`, which is then linked with the `M7` program executable and the `unittest` GoogleTest executable.

For production calculations on HPC clusters, the native installation procedure is recommended.

The native procedure is generally straightforward if the user simply wishes to try out the program on small-scale 
calculations or investigate some potential algorithmic improvement. Additionally however, a containerized solution is
detailed in the "Docker image" section below as "run anywhere" alternative suited to this exploratory use case.

## Native build
M7 requires:
* C++11 compiler implementing MPI
* CMake installation at least as recent as version 3.10
* LAPACK implementation

First, retrieve the program code from this repository:
```bash
git clone https://github.com/robert-anderson/M7.git
```
Then create and `cd` to a directory `<BIN_PATH>` in which CMake can build its targets and run `cmake` with MPI compilers
```bash
mkdir -p <BIN_PATH>; cd <BIN_PATH>
CC=mpicc CXX=mpicxx FC=mpif90 cmake <SRC_PATH> <OPTS>
```
Where `<SRC_PATH>` should be replaced with the path of the cloned repository, and `<OPTS>` with all desired CMake compile time variables. 
Finally, build all targets in parallel with `make -j`. The program binary will then be found at `bin/M7`, and that of the unit tests at `test/unittest`.

N.B. If using the Intel Parallel Studio compilers, make sure the compilevars.sh has been sourced to correctly configure the library path environment variables prior to invoking CMake.
Otherwise, HDF5 may fail to find Intel symbols such as `_intel_fast_memcpy` at link time.

## Docker image
This may be a useful alternative to the native procedure for users on a variety of different operating systems and architectures.
A Docker installation is required for this procedure (https://docker.com).
At the command line, pull the latest version of the image and tag it with a locally-resolved name
```bash
docker pull ghcr.io/robert-anderson/m7_tester:latest
docker tag ghcr.io/robert-anderson/m7_tester:latest m7_tester
```
Now a terminal session can be accessed within a container instantiated from the downloaded image:
```bash
docker run -it m7_tester bash
```
The session will begin as root within the container. Since some MPI implementations warn against running as root, M7 
should be built under a non-root account. The Docker image includes one called `user`:
```bash
su user
cd ~
```
From this point, one can follow the native build procedure, or run the included scripts to automate the process. For example:
```
~/build.sh Release fermion no 4
```
will checkout the latest version of this repo before building the `Release` version of `M7` and `unittest`, with fermion `MBF_TYPE`, no complex arithmetic (real integrals and walkers), and parallel compilation over 4 threads.  

# Running the test suite
Testing is split into:
* Unit tests, which use the GoogleTest framework to verify individual classes and functions in a special test-specific binary called `unittest`.
* System tests, which runs the `M7` program binary and performs two kinds of checks:
  * Comparative: checks current output against that of a trusted reference run 
  * Static: checks current output against some externally-defined correct values

Unit tests can be run in full by simply executing `unittest` for verification of serial (1 MPI rank) runs, or 
`mpirun -n <NRANK> unittest` to verify for any number of MPI ranks `<NRANK>`. If the user is only interested in a subset
of tests, this can be indicated in the `--gtest-filter` option, e.g.
```bash
unittest --gtest_filter=FrmOnvField.*
```
will only run the unit tests within the `FrmOnvField` section which tests all features of M7's Slater determinant encoding class.

System tests are implemented in Python.
A system test can have the following outcomes:
* PASS: all comparative and static checks passed for the current run
* FAIL: at least one comparative or static check failed for the current run
* SKIP: the binary was not built with compile time options compatible with this test definition, or the test script contains a python error

A single test can be run in isolation just by executing its test script e.g.:
```bash
cd <SRC_PATH>/system_test
export PYTHONPATH=PYTHONPATH:$(pwd)
python ./frm/real/abinitio/energy/uhf/parallel/test.py <BIN_PATH>/bin/M7
```
will run the system test which verifies the correctness of estimating a ground state energy with integrals derived from 
an unrestricted Hartree-Fock calculation with an MPI-parallel M7 run.
This procedure can be extended to running several system tests in parallel using the `run.py` script, which ensures that
the MPI slots on the testing machine are never oversubscribed. e.g.
```bash
python run.py test <BIN_PATH>/bin/M7 mpirun -p $(find . -name test.py)
```
will find all system test scripts and execute them, performing both static and comparative (if reference is defined) checks
and printing out the `PASS`/`FAIL`/`SKIP` outcome of each job as they terminate. 
The `find` call may of course be replaced by a list of script paths, or any other command the user may decide to use in the
selection of a subset of tests to run. Additionally, with the command
```bash
python run.py redef <BIN_PATH>/bin/M7 mpirun -p <LIST_OF_SCRIPT_PATHS>
```
the user signals the intent to redefine the references for every system test in the supplied list. In this case, only 
the static checks are performed, and if every test in the set passes, the output artefacts defining the reference runs are
updated.

Since the system tests are just automated executions of the program, the defintions within the `system_test`
directory are a source of example `config.yaml` scripts on which users may wish to base their own input files.

All unit and system tests are executed by a GitHub Action on commits and PRs to `dev` branch on this repo.


# Usage
The program takes all run time variables by parsing a YAML configuration file. To run a calculation specified by `config.yaml`, simply pass the path to `config.yaml` as the only command line argument to the M7 binary. 
Otherwise, if no argument is specified, the complete specification for the configuration document is written to standard output.
Piping this output to `less -RS` may aid the user in searching for particular options.