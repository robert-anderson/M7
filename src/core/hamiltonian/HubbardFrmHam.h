//
// Created by anderson on 12/8/21.
//

#ifndef M7_HUBBARDFRMHAM_H
#define M7_HUBBARDFRMHAM_H

#include <src/core/util/Foreach.h>
#include <src/core/nd/NdFormatD.h>
#include "FrmHam.h"
#include "src/core/linalg/Sparse.h"

struct HubbardFrmHam : FrmHam {
    /**
     * multidimensional format describing the orthogonally-coordinated lattice sites
     */
    const NdFormatD m_format;
    /**
     * boundary conditions for each lattice dimension (0 for OBC, 1 for PBC, -1 for APBC)
     */
    const std::vector<int> m_bcs;
    /**
     * on-site repulsion scalar in units of the hopping
     */
    const defs::ham_t m_u;
    /**
     * sparse map enabling lookup of all coordinated sites given a row site
     */
    sparse::Matrix<int> m_t_mat_sparse;
    /**
     * dense map of coordinated sites allowing lookup of the 1-body H matrix element given the flat indices of two site
     * index vectors
     */
    dense::Matrix<int> m_t_mat_dense;
    /**
     * to help with excitation generation.
     */
    size_t m_unique_nconn_product;
    /**
     * whether the model meets the sign problem-free conditions
     */
    const bool m_spf;
private:

    size_t get_coord_index(const defs::inds &site_inds, size_t idim, size_t value) const;

    std::pair<size_t, int> get_coordination(const defs::inds &site_inds, size_t idim, bool inc) const;

    static size_t nsite(const defs::inds& site_shape);

    /**
     * determine whether this hubbard hamiltonian is sign-problematic
     * @return
     *  true if the model is sign problem-free
     */
    bool sign_problem() const;

public:

    HubbardFrmHam(const defs::inds& site_shape, const std::vector<int>& bcs,
                       defs::ham_t u, int ms2_restrict, int charge);

    HubbardFrmHam(const fciqmc_config::FermionHamiltonian& opts);

    defs::ham_t get_coeff_1100(const size_t &i, const size_t &j) const override;

    defs::ham_t get_coeff_2200(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const override;


    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    void log_data() const override;
};


#endif //M7_HUBBARDFRMHAM_H
