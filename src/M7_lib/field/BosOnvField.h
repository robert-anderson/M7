//
// Created by Robert J. Anderson on 07/04/2021.
//

#ifndef M7_BOSONVFIELD_H
#define M7_BOSONVFIELD_H

#include <M7_lib/basis/BasisData.h>
#include <M7_lib/caches/DecodedBosOnv.h>

#include "NumberField.h"

/*
 * forward declarations to support occupation factor methods
 */
struct BosOnvConnection;
class BosOps;


struct BosOnvField : NdNumberField<defs::bos_occ_t, 1> {
    typedef NdNumberField<defs::bos_occ_t, 1> base_t;
    /**
     * single-particle basis data is all the system-related data that this class needs to retain
     */
    const sys::bos::Basis m_basis;
    /**
     * a refreshable cache of useful representations for excitation generation and enumeration
     */
    mutable decoded_mbf::BosOnv m_decoded;

    BosOnvField(Row *row, const sys::bos::Basis &basis, std::string name = "");

    /*
     * all ONVs implement the following ctor
     */
    BosOnvField(Row *row, const sys::Basis &basis, std::string name = "");

    /*
     * this particular MBF only needs the basis, but future MBF types might need the full sector information, and so
     * a common interface is realised by implementing a ctor of the following form in all MBFs
     */
    BosOnvField(Row *row, const sys::bos::Sector &sector, std::string name = "");
    BosOnvField(Row *row, const sys::Sector &sector, std::string name = "");

    BosOnvField(const BosOnvField &other);

    BosOnvField &operator=(const BosOnvField &other) {
        base_t::operator=(other);
        return *this;
    }

    BosOnvField &operator=(const defs::uintv_t &inds);

    /**
     * set based on the (repeated) boson operator indices
     * @param iops
     *  operator indices (one for each boson in the state)
     */
    void set_ops(const defs::uintv_t &iops);

    size_t nboson() const;

    /**
     * compute the "occupation factor" required to keep the boson ONV basis orthonormal.
     * b |n>  = sqrt(n) |n-1>
     * b+ |n> = sqrt(n+1) |n+1>
     *
     * b^m |n>  = sqrt(n(n-1)...(n-m+1)) |n-m>
     * b+^m |n>  = sqrt((n+1)(n+2)...(n+m)) |n+m>
     * @param src
     * @return
     */
    size_t occ_fac_square(const BosOnvConnection &conn) const;

    double occ_fac(const BosOnvConnection &conn) const;

private:
    /**
     * @param occ
     *  occupation of the mode
     * @param nop
     *  power of the annihilation operator acting on the mode
     * @return
     *  square of the occupation factor associated with the annihilation
     */
    static size_t occ_fac_square_ann(size_t occ, size_t nop);

    /**
     * @param occ
     *  occupation of the mode
     * @param nop
     *  power of the creation operator acting on the mode
     * @return
     *  square of the occupation factor associated with the creation
     */
    static size_t occ_fac_square_cre(size_t occ, size_t nop);

    /**
     * @param occ
     *  occupation of the mode
     * @param nop
     *  power of the number operator acting on the mode
     * @return
     *  square of the occupation factor associated with the number operator
     */
    static size_t occ_fac_square_com(size_t occ, size_t nop);

public:

    /**
     * if the matrix element in question is actually contracted, with common indices among the creation and annihilation
     * operators, the occupation factors change.
     * the ONV can be rearranged to group operators in mode-wise normal order. for a given mode, common indices result
     * in the following factor:
     * b+b |n> = n |n>
     * b+^2b^2 |n> = n(n-1) |n>
     * b+^mb^m |n> = n(n-1)...(n-m+1) |n>
     * @param src
     * @param com
     * @return
     */
    size_t occ_fac_square(const BosOnvConnection &conn, const BosOps &com) const;

    double occ_fac(const BosOnvConnection &conn, const BosOps &com) const;

    /**
     * same as calling occ_fac_square on an empty BosOnvConnection (but avoids the need to allocate the connection)
     * @param com
     *  indices representing pairs of creation and annihilation ops
     * @return
     *  diagonal occupation factor
     */
    size_t occ_fac_square(const BosOps &com) const;

    double occ_fac(const BosOps &com) const;
};


#endif //M7_BOSONVFIELD_H
