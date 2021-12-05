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
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j M7
```
The binary will then appear at src/M7. `Release` can be replaced with `Debug` above to build the debug binary.

Parallel HDF5 is a non-prerequisite dependency which is downloaded, compiled, and linked as a static library during the build process.
If using the Intel C++ compiler, make sure the compilevars.sh has been sourced to correctly configure the library path environment variables prior to invoking CMake.
Otherwise, HDF5 may fail to find Intel symbols such as `_intel_fast_memcpy` at link time.

Also, the HDF5 build script may quit, reporting that `MPI_C_FOUND` is not set, meaning that CMake's `FIND_PACKAGE` feature could not determine the path of the MPI library.
Ensuring that the MPI library root is included in the `PATH` environment variable should be sufficient to overcome this.
This issue may arise in IDEs where the MPI compiler executable has been configured as the default compiler, which allows the M7 CMakeLists.txt file to find the MPI path, but fails for the HDF5 external project.
In such cases, set the PATH variable within the IDE's toolchain dialogue to contain the MPI library root.


# Usage
The program takes all run time variables through the command line interface.
The command line options available can be listed by executing the debug or release binaries with the -h or --help option.

Complex configurations can become cumbersome to manage on the command line directly, so particular run time configurations can easily be split over multiple lines in a shell script.

If the user chooses, the python/ directory contains a simple and convenient wrapper which encapsulates all options, and can handle the execution of distributed instances of M7.


# General remarks for contributing developers

## Symbol conventions

Unless for temporary and simple loop dummy indices, all symbols are to be descriptive.
Many mature IDEs are available for C++11 which have symbol completion and other productivity enhancing features, and it is assumed (though not strictly required) that contributors will be make use one of these programs.

CamelCase is the convention for classes, and snake_case for methods and members.
Struct/class member symbols are always prefixed with `m_`, and constant expressions with `c_`.
Type definitions are suffixed with `_t`, but a single uppercase character is commonly used for template arguments.

Global data of any kind is strongly discouraged, but in the rare circumstances that global symbols are warranted (e.g. coupling to an external dependency), the prefix `g_` is used.

The private/protected members + getters and setters style of member access leads to clutter and is discouraged, and const correctness is instead emphasised (i.e. a read-only member should be public and const, constructed by a "make" method if required).
Declaring methods within private/protected blocks on the other hand, increases code readability by clarifying intent in most cases, and is certainly encouraged.

Validation macros `ASSERT` and `REQUIRE` are important in ensuring run-time correctness.
All such macros containing `ASSERT` are only compiled in builds without the preprocessor symbol `NDEBUG` defined.
Thus, these macros should be instated liberally, for example in bounds checking situations so that overheads associated with such validations can be omitted in the `release` build.

The `REQUIRE` macros are compiled in every build target, and behave like uncaught exceptions.
They should be used only in rarely-called routines, e.g. checking for existence of file in setup of I/O-performing object.

Compile-time correctness should be enforced with the compiler's builtin `static_assert` method.

With regard to parallelization, there are three kinds of method:
* Independent: makes no internal reference to MPI symbols.
* Collective: makes reference to MPI symbols and must be called by all ranks.
* Selective: is called by a collective method on only a subset of ranks.
If a method is MPI selective, its symbol is to be suffixed with an underscore.
  <!---TODO parallelization in concert with validation macros -->
  
Although this project makes heavy use of template metaprogramming, wherever possible templating should be avoided in situations where function overloading could obtain the same result.
Untemplated method bodies in untemplated classes should be specified in the implementation file, except for any that are very frequently called - the compiler cannot inline methods when only the prototype is available in the header. 

Templated classes require the explicit declaration of specializations in order for the method definitions to be relocated from the `.h` header file to the `.cpp` implementation file, and this is not a practice followed in this project. 
Hence, all templated code will remain in the header files.

Developers should be aware that while it is a tidy approach, runtime polymorphism (virtual-override idiom) is not a zero-cost abstraction, and therefore it is not advised for frequently invoked methods.
The low-level "datasystem" modules have avoided virtual methods, which should be the most frequently called, so in the development of future functionality, the use of this polymorphism is likely beneficial in terms of code readability and will not degrade performance noticeably.

C preprocessing directive blocks should never appear outside of `src/defs.h` and the few other places where macros are defined.
However, when introducing a new piece of optional functionality, it is encouraged that a corresponding `ENABLE_` macro be introduced, along with a `constexpr` boolean variable which is to be used in other parts of the code to effect compile-time branching.
Otherwise, new modules cause unnecessary run-time branching, potentially having an adverse effect on performance.