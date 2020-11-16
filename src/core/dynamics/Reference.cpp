//
// Created by rja on 03/07/2020.
//

#include "Reference.h"
#include "Wavefunction.h"

Reference::Reference(Wavefunction &wf, const Hamiltonian& ham, views::Onv &onv,
                     const Options &m_opts) :
        elements::Onv(onv), m_wf(wf), m_ham(ham), m_aconn(onv),
        m_redefinition_thresh(m_opts.reference_redefinition_thresh){

    m_irank = m_wf.m_ra.get_rank(onv);
    /*
     * check for onv in walker table first
     */
    auto irow = wf.m_walkers[onv];
    m_irow = ~0ul;
    if (irow) {
        m_irow = *irow;
        ASSERT(m_irank==mpi::irank());
    }
    else if (mpi::i_am(m_irank)) {
        m_irow = m_wf.create_walker(onv, m_opts.nwalker_initial, ham.get_energy(onv), true, true);
    }
    change(m_irow, m_irank);
}