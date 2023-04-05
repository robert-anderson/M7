//
// Created by rja on 03/11/22.
//

#ifndef M7_FOCKRDM4_H
#define M7_FOCKRDM4_H

#include "Rdm.h"

struct FockMatrix : dense::SquareMatrix<ham_t> {
    FockMatrix(uint_t nsite, const str_t& fname): dense::SquareMatrix<ham_t>(nsite) {
        hdf5::FileReader reader(fname);
        REQUIRE_TRUE(reader.child_exists("ACT_FOCK_INDEX"), "invalid fock matrix file contents");
        REQUIRE_TRUE(reader.child_exists("ACT_FOCK_VALUES"), "invalid fock matrix file contents");
        const auto inds = hdf5::DatasetLoader::load_vector<int64_t>(reader, "ACT_FOCK_INDEX");
        const auto values = hdf5::DatasetLoader::load_vector<double>(reader, "ACT_FOCK_VALUES");
        REQUIRE_EQ(inds.size(), values.size()*2, "incorrect number of indices");
        for (uint_t i = 0ul; i<values.size(); ++i) {
            // offset by one to account for the Fortran indices in the HDF5 format
            const auto irow = inds[2*i]-1;
            const auto icol = inds[2*i+1]-1;
            REQUIRE_GE(irow, 0, "Fock matrix row index OOB (HDF5 indices dataset must count from 1)");
            REQUIRE_GE(icol, 0, "Fock matrix col index OOB (HDF5 indices dataset must count from 1)");
            (*this)(irow, icol) = values[i];
        }
        symmetrize();
    }
};

class FockRdm4 : public ContractedRdm {
protected:
    /**
     * if the Fock matrix has no non-zero diagonal elements, then F4RDM takes no contributions from walker death or
     * block averaging
     */
    const bool m_nonzero_diagonal;
public:
    FockRdm4(const conf::Rdms &opts, OpSig max_contrib_exsig, sys::Sector sector,
             uint_t nvalue, bool nonzero_diagonal);
};

class NonDiagFockRdm4 : public FockRdm4 {
    /**
     * Generalized Fock in the active space
     */
    const dense::SquareMatrix<ham_t> m_fock;

public:

    NonDiagFockRdm4(const conf::Rdms &opts, const FockMatrix& fock, sys::Sector sector, uint_t nvalue);

    /**
     * override the default method to implement on-the-fly contraction
     */
    void frm_make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn,
                           const FrmOps& com, wf_t contrib) override;
};

class DiagFockRdm4 : public FockRdm4 {
    /**
     * Diagonal elements of the generalized Fock in the active space
     */
    const dense::Vector<ham_t> m_fock;

public:

    DiagFockRdm4(const conf::Rdms &opts, const FockMatrix& fock, sys::Sector sector, uint_t nvalue);

    /**
     * override the default method to implement on-the-fly, diagonal-only contraction
     */
    void frm_make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn,
                           const FrmOps& com, wf_t contrib) override;
};

#endif //M7_FOCKRDM4_H
