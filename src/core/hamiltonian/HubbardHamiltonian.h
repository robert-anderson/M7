//
// Created by anderson on 12/8/21.
//

#ifndef M7_HUBBARDHAMILTONIAN_H
#define M7_HUBBARDHAMILTONIAN_H

#include <src/core/util/Foreach.h>
#include <src/core/nd/NdFormatD.h>
#include "FermionHamiltonian.h"
#include "src/core/linalg/Sparse.h"

struct HubbardHamiltonian : FermionHamiltonian {
    const NdFormatD m_format;
    const std::vector<int> m_bcs;
    const defs::ham_t m_u;
    sparse::Matrix<int> m_t_mat_sparse;
    Matrix<int> m_t_mat_dense;
private:

    size_t get_coord_index(const defs::inds &site_inds, size_t idim, size_t value) const;

    std::pair<size_t, int> get_coordination(const defs::inds &site_inds, size_t idim, bool inc) const;

    static size_t nsite(const defs::inds& site_shape);

public:

    HubbardHamiltonian(const defs::inds& site_shape, std::vector<int> bcs, defs::ham_t u, int ms2_restrict, int charge);

    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    void log_data() const override;
};


#endif //M7_HUBBARDHAMILTONIAN_H
