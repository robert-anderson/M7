//
// Created by rja on 03/11/22.
//

#ifndef M7_FOCKRDM4_H
#define M7_FOCKRDM4_H

#include "Rdm.h"

struct FockMatrix : dense::SquareMatrix<ham_t> {
    FockMatrix(uint_t nsite, const str_t& fname): dense::SquareMatrix<ham_t>(nsite) {
        logging::info("Loading generalized Fock matrix for accumulation of its contraction with the 4RDM");
        hdf5::FileReader reader(fname);
        REQUIRE_TRUE(reader.child_exists("ACT_FOCK_INDEX"), "invalid fock matrix file contents");
        REQUIRE_TRUE(reader.child_exists("ACT_FOCK_VALUES"), "invalid fock matrix file contents");
        v_t<int64_t> inds;
        reader.read_data("ACT_FOCK_INDEX", inds);
        v_t<double> values;
        reader.read_data("ACT_FOCK_VALUES", values);
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

class FockRdm4 : public Rdm {
    /**
     * Generalized Fock in the active space
     */
    const dense::SquareMatrix<ham_t> m_fock;

public:

    FockRdm4(const conf::Rdms &opts, const FockMatrix& fock, sys::Sector sector, uint_t nvalue);

    /**
     * override the default method to implement on-the-fly contraction
     */
    void frm_make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn,
                           const FrmOps& com, const wf_t& contrib) override;
};

//class DiagFockRdm4 : public Rdm {
//    /**
//     * Diagonal elements of the generalized Fock in the active space
//     */
//    v_t<ham_t> m_fock;
//
//public:
//    /**
//     * assume diagonal Fock matrix until a non-diagonal value is read-in
//     */
//    bool m_diagonal;
//
//    FockRdm4(const conf::Rdms& opts, sys::Size basis_size, uint_t nelec, uint_t nvalue);
//
//    /**
//     * override the default method to implement on-the-fly contraction
//     */
//    void frm_make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn,
//                           const FrmOps& com, const wf_t& contrib) override;
//};

#endif //M7_FOCKRDM4_H
