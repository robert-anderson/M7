//
// Created by Robert J. Anderson on 07/04/2021.
//

#ifndef M7_BOSONVFIELD_H
#define M7_BOSONVFIELD_H

#include <M7_lib/basis/BasisData.h>
#include <M7_lib/caches/DecodedBosOnv.h>

#include "NumberField.h"
#include "M7_lib/foreach/SetByteForeach.h"

/*
 * forward declarations to support occupation factor methods
 */
struct BosOnvConnection;
class BosOps;


struct BosOnvField : NdNumberField<bos_occ_t, 1> {
    typedef NdNumberField<bos_occ_t, 1> base_t;
    /**
     * single-particle basis data is all the system-related data that this class needs to retain
     */
    const sys::bos::Basis m_basis;
    /**
     * a refreshable cache of useful representations for excitation generation and enumeration
     */
    mutable decoded_mbf::BosOnv m_decoded;

    BosOnvField(Row *row, const sys::bos::Basis &basis, str_t name = "");

    /*
     * all ONVs implement the following ctor
     */
    BosOnvField(Row *row, const sys::Basis &basis, str_t name = "");

    /*
     * this particular MBF only needs the basis, but future MBF types might need the full sector information, and so
     * a common interface is realised by implementing a ctor of the following form in all MBFs
     */
    BosOnvField(Row *row, const sys::bos::Sector &sector, str_t name = "");
    BosOnvField(Row *row, const sys::Sector &sector, str_t name = "");

    BosOnvField(const BosOnvField &other);

    BosOnvField &operator=(const BosOnvField &other) {
        base_t::operator=(other);
        return *this;
    }

    BosOnvField &operator=(const uintv_t &inds);

    /**
     * set based on the (repeated) boson operator indices
     * @param iops
     *  operator indices (one for each boson in the state)
     */
    void set_ops(const uintv_t &iops);

    uint_t nboson() const;

    template<typename fn_t>
    void foreach_setmode(const fn_t& fn) const {
        static_assert(std::is_same<bos_occ_t, buf_t>::value,
                "this approach uses aliasing, so the boson occupation must have the same type as the buffer");
        auto get_work_fn = [this](uint_t idataword){
            return reinterpret_cast<uint_t*>(begin())[idataword];
        };
        setbyte_foreach::single<uint_t>(integer::divceil(m_size, sizeof(uint_t)), fn, get_work_fn);
    }

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
    uint_t occ_fac_square(const BosOnvConnection &conn) const;

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
    static uint_t occ_fac_square_ann(uint_t occ, uint_t nop);

    /**
     * @param occ
     *  occupation of the mode
     * @param nop
     *  power of the creation operator acting on the mode
     * @return
     *  square of the occupation factor associated with the creation
     */
    static uint_t occ_fac_square_cre(uint_t occ, uint_t nop);

    /**
     * @param occ
     *  occupation of the mode
     * @param nop
     *  power of the number operator acting on the mode
     * @return
     *  square of the occupation factor associated with the number operator
     */
    static uint_t occ_fac_square_com(uint_t occ, uint_t nop);

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
    uint_t occ_fac_square(const BosOnvConnection &conn, const BosOps &com) const;

    double occ_fac(const BosOnvConnection &conn, const BosOps &com) const;

    /**
     * same as calling occ_fac_square on an empty BosOnvConnection (but avoids the need to allocate the connection)
     * @param com
     *  indices representing pairs of creation and annihilation ops
     * @return
     *  diagonal occupation factor
     */
    uint_t occ_fac_square(const BosOps &com) const;

    double occ_fac(const BosOps &com) const;
};


#endif //M7_BOSONVFIELD_H
