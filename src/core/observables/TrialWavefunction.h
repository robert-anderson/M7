/**
 * @file
 * @author Robert John Anderson <robert.anderson@kcl.ac.uk>
 *
 * @section LICENSE
 *
 * @section DESCRIPTION
 * Trial wavefunction is the standard procedure for extracting Hamiltonian
 * Rayleigh quotient numerator and denominator time series from an FCIQMC
 * walker distribution. Denoting the trial wavefunction \f$\ketTrial\f$, the
 * numerator is \f$\braTrialn\Hop\ketPsin\f$, and the denominator is the
 * overlap \f$\braTrialn\Psin\rangle\f$.
 *
 * The single determinant case is handled in the Reference class, since the
 * connection of a given occupied determinant in the walker list to the
 * reference is trivially obtainable from the excitation level.
 *
 * The multideterminant case is introduced in order to reduce the random error
 * in the numerator and denominator time series and consequently arrive at a more
 * precise statistical energy estimation. However, in this case, ascertaining
 * whether a given occupied determinant is connected to the trial wavefunction
 * subspace is more involved, and efficient evaluation of the correct numerator
 * contribution requires precomputation. Dropping the cycle index \f$n\f$, and
 * denoting the trial and connected determinant subspaces as \f$\mathcal{T}\f$
 * and \f$\mathcal{C}\f$ we can write the numerator expression in the many-body basis
 * \f[
 * \braTrial\Hop\ketPsi = \sum_{\bfi \in \mathrm \mathcal{T}}
 *      \CiTrial^* \braDi \Hop \sum_{\bfj \in \mathrm \mathcal{C}} \Cj \Dj
 * \f]
 * Which can be reorganised as a simple scalar product between the current
 * walker populations in the connection space with a vector containing the
 * projection of the hamiltonian onto the static trial wavefunction.
 * \f[
 * \braTrial\Hop\ketPsi = \sum_{\bfj \in \mathrm \mathcal{C}}
 *      \left( \sum_{\bfi \in \mathrm \mathcal{T}} \CiTrial^* \Hij \right) \cdot \Cj
 * \f]
 * The denominator is evaluated simply by summing together the conjugated
 * product of trial wavefunction amplitudes and walker populations.
 * * \f[
 * \braTrialn\Psin\rangle = \sum_{\bfi \in \mathrm \mathcal{T}}
 *      \CiTrial^* \braDi \Hop \sum_{\bfj \in \mathrm \mathcal{C}} \Cj \Dj
 * \f]
 *
 * The other feature which distinguishes the general case of trial wavefunction
 * projection from the Reference case, is the need to solve the correlation problem
 * within the subspace. In this class, the local trial subspace is set up, then
 * the build method gathers the full subspace on every process. The root process
 * then diagonalises the DenseHamiltonian defined by the full DeterminantList.
 *
 * The appropriate eigenvectors are then broadcast to all ranks, and the local
 * connected space is constructed. The connected space is a MappedDeterminantList,
 * because this allows the connected space to be probed in constant time.
 *
 */

#ifndef M7_TRIALWAVEFUNCTION_H
#define M7_TRIALWAVEFUNCTION_H


#include "src/core/dynamics/WalkerList.h"
#include "src/core/basis/DeterminantList.h"
#include "src/core/basis/MappedDeterminantList.h"

class TrialWavefunction {

    WalkerList *m_walker_list = nullptr;
    /*
    DeterminantList m_local_subspace;
    //DeterminantList m_full_subspace;
    MappedDeterminantList m_connected_space;

    std::vector<defs::wf_t> m_local_weights;
    */


public:
    TrialWavefunction(WalkerList *walker_list) : m_walker_list(walker_list){

    }

    void add_determinant(const size_t& irow){

    }

};


#endif //M7_TRIALWAVEFUNCTION_H
