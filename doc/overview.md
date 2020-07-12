# Overview

M7 (*Many-body Stochastic Expectation Value Estimation Networks*) is a stochastic quantum chemistry software package which primarily implements the FCIQMC method \cite doi:10.1063/1.3193710.
The purpose of this overview page is to help users and developers understand where each module fits into the program at large.
Modules will be introduced in a section on this page with a coarse-grained exposition of their relationships with other parts of the program.
Each section will provide at least one hyperlink to the C++ class in which each piece of referenced functionality is implemented.
The linked-to Doxygen pages will display the detailed information about the class, its members and methods, and thorough instructions for its utilization.
Thus the documentation is almost completely code-driven, with the responsibility for maintaining the docs for a given module lying with the developer(s) who wrote the corresponding header files.
The intent of this documentation website is to provide a space in which M7 users and developers can learn from and contribute to a body of knowledge, advice, and best practices accumulated over the entire lifetime of the project.

## FCIQMC Algorithm

The FCIQMC method is founded on the assumption that stochastic application of the projector method recursion relation
\f[
    \label{eq:master}
    \ketPsinext = (1-\tau\Hop)\ketPsin
\f]
with a suitably small *timestep* \f$\tau\f$, 
where the \ref Wavefunction is a linear superposition of Slater determinants
\f[
    \label{eq:psidef}
    \ketPsi \equiv \sum_\bfi \Ci \ketDi
\f]
can yield accurate eigenvalues and other properties for configuration-interation problems far beyond the scope of exact diagonalisation methods.

This delay in the onset of the "curse of dimensionality" faced by exponentially-scaling correlation methods is achieved by discretising the \f$\Ci\f$ coefficients as a population of *walkers*.
Thus, at any given iteration of \f$\eqref{eq:master}\f$, the number of determinants with any walker occupation is usually a shrinking minority of the total dimension of the space (provided that a single particle basis conducive to this sparsity e.g. the canonical Hartree-Fock orbitals is chosen).
By this approach, the memory scaling of the FCI problem can be brought within practicable limits, provided that an efficient system of data structures can be devised to store and update a sparse, parallelised representation of \f$\ketPsi\f$

Computational tractability is also aided by the walker discretisation, since the accumulation of walkers on a given determinant indicates the number of attempts the algorithm should make to convey new walkers from that source.
The processes by which the \f$\Ci\f$ coefficients are updated are called *spawning* (off-diagonal elements), and *death* (off-diagonal elements)
\f[
    \Cinext = \Cin - \tau \sum_{\bfj\neq\bfi} \Hij \Cjn - \tau(\Hii - S) \Cin
\f]
where \f$S\f$ is an approximation to the exact ground-state eigenvalue called the *diagonal shift*.
In the expression of the FCIQMC algorithm, it is convenient to use the following notation
\f[
    \Cinext \equiv \Cin + \dspawnn + \ddeathn
\f]

Structurally, the FCIQMC code is contained within a \ref FciqmcCalculation object, which serves as a high-level class tying together the more intricate objects in the implementation.
One such crucial object is the \ref Wavefunction, which in turn contains a system of three central \ref List objects: namely a \ref WalkerList --- which records the state of the discretised \f$\ketPsin\f$ at an iteration \f$n\f$, and two \ref SpawnList objects --- which are responsible for the MPI-based communication of stochastically-generated off-diagonal Monte Carlo moves.

### Spawning and Death
 walker list            | send buffer    | receive buffer |
------------------------|----------------|----------------|
 \f$\Cin\f$             | \f$0\f$        | \f$0\f$        |
The walker list is sequentially accessed, with each *source* \ref Determinant \f$\Dj\f$ selecting zero or more *destination* determinants through an \ref ExcitationGenerator.
The instantaneous walker weight on \f$\Dj\f$ is conveyed as the ratio of \f$\Cjn\f$ to the number of attempts made to generate a connected \f$\Di\f$.
The death contribution is taken into account in the same loop after all spawns have been generated. 
 walker list            | send buffer    | receive buffer |
------------------------|----------------|----------------|
 \f$\Cin+\ddeathn\f$    | \f$\dspawnn\f$ | \f$0\f$        |
at the end of this loop over occupied determinants, the walker list reflects only the update to the approximate wavefunction due to the diagonal part of the shifted Hamiltonian.

### Communication
The generated spawns have been enumerated in a segmented \ref List, with the segment corresponding to the id of the MPI rank due to receive the spawned contribution. This rank is decided on a determinant block-by determinant block basis through a \ref RankAllocator, which can dynamically reallocate determinant blocks between MPI ranks an effort to eliminate performance stifling load imbalance.

Sent spawns are scattered by and gathered to every process via the MPI_Alltoallv subroutine.

 walker list            | send buffer    | receive buffer |
------------------------|----------------|----------------|
 \f$\Cin+\ddeathn\f$    | \f$0\f$        | \f$\dspawnn\f$ |
### Annihilation
This step is a sequential access on the incoming spawned contributions, resulting in random access update of the walker list by hash table lookup.

 walker list            | send buffer    | receive buffer |
------------------------|----------------|----------------|
 \f$\Cinext\f$          | \f$0\f$        | \f$0\f$        |


## Hybrid Parallelized Variables
Suppose we have an integer variable `n` which is accumulated in a hybrid parallel scheme which preserves the value of `n` local to the MPI rank.
That is to say that the final value of `n` is an `MPI_SUM` reduction of a sum reduction over partial values on each thread.
This leads to an ugly proliferation of symbols e.g.
```cpp
size_t n = 0ul; // the final MPI-reduced integer
size_t n_local = 0ul; // the local thread-reduced integer
#pragma omp parallel default(none)
{
    size_t n_thread = do_something(); // temporary variable private to the thread
#pragma omp atomic update
    n_local+=n_thread;
}
MPI_Allreduce(&n_local, &n, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_WORLD)
```
We are in need of a solution to this problem in order to keep the code maintainable.
In this example situation, the \ref Hybrid class would be used to encapsulate all relevant quantities
```cpp
Hybrid<size_t> n;
#pragma omp parallel default(none)
{
    n.thread() = do_something(); // temporary variable private to the thread
}
n.sum(); // returns result of hybrid reduction
```
In situations where there is no need for thread-reduction, the base class \ref Reducable is an adequate solution.


<!--
In the open-source domain, there are two notable FCIQMC implementations which predate this project, namely \link https://github.com/ghb24/NECI_STABLE NECI\endlink and \link https://github.com/hande-qmc HANDE\endlink.
Of these two projects, the most mature and fully-featured is NECI. It is also the mose productive in terms of research output, and is the code with which the Booth group at King's College London is most familiar.
For these reasons, the following mission statement of the present project shall be stated in terms that reference certain perceived shortcomings of the NECI implementation, while acknowledging the vast effort invested into making NECI such a capable and efficient implementation.

The development of M7 aims to:
* **Set a high standard for code quality**
* **Encapsulate the modern FCIQMC algorithm**
* **Maintain a well-structured codebase**
-->
