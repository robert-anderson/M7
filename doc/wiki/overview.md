# Wiki

M7 (*Many-body Stochastic Expectation Value Estimation Networks*) is a stochastic quantum chemistry software package which primarily implements the FCIQMC method \cite doi:10.1063/1.3193710.
The purpose of this wiki is to provide a space in which M7 users and developers can learn from and contribute to a body of knowledge, advice, and best practices accumulated over the entire lifetime of this project.

The FCIQMC method is founded on the assumption that stochastic application of the projector method recursion relation
\f[
    \label{eq:master}
    \ketPsinext = (1-\tau\Hop)\ketPsin
\f]
where the wavefunction is a linear superposition of Slater determinants
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


<!--
In the open-source domain, there are two notable FCIQMC implementations which predate this project, namely \link https://github.com/ghb24/NECI_STABLE NECI\endlink and \link https://github.com/hande-qmc HANDE\endlink.
Of these two projects, the most mature and fully-featured is NECI. It is also the mose productive in terms of research output, and is the code with which the Booth group at King's College London is most familiar.
For these reasons, the following mission statement of the present project shall be stated in terms that reference certain perceived shortcomings of the NECI implementation, while acknowledging the vast effort invested into making NECI such a capable and efficient implementation.

The development of M7 aims to:
* **Set a high standard for code quality**
* **Encapsulate the modern FCIQMC algorithm**
* **Maintain a well-structured codebase**
-->
