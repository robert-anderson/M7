//
// Created by rja on 26/10/22.
//

#ifndef M7_OPSIG_H
#define M7_OPSIG_H

#include "M7_lib/util/Integer.h"

struct OpSig {
    /**
     * number of bits in the signature representing each number of fermion SQ operators
     */
    static constexpr uint_t c_nbit_nop_frm = 3;
    /**
     * number of bits in the signature representing each number of boson SQ operators
     */
    static constexpr uint_t c_nbit_nop_bos = 2;
    /**
     * mask and max value for extraction of a number of fermion SQ operators
     */
    static constexpr uint_t c_nop_mask_frm = (1 << c_nbit_nop_frm) - 1;
    /**
     * mask and max value for extraction of a number of boson SQ operators
     */
    static constexpr uint_t c_nop_mask_bos = (1 << c_nbit_nop_bos) - 1;
    /**
     * total number of distinct excitation signatures that can be stored
     */
    static constexpr uint_t c_ndistinct = (1 << (2 * c_nbit_nop_frm + 2 * c_nbit_nop_bos));

private:
    const uint_t m_i;

    /**
     * compactly expresses an arbitrary SQ operator product as a single integer given some compile-time constant numbers
     * of bits for each element. e.g. if c_nbit_nop_frm = 3 and c_nbit_nop_bos = 1, then 2x3+2x1 = 8 bits are
     * required to store a connection excitation level as an exsig (excitation signature) with upto 7 fermion creation
     * operators, 7 fermion annihilation operators and 1 each of boson creation and annihilation operators, this limit
     * should be sufficient for all foreseeable applications, but these bit segment lengths are not hardcoded.
     * @param nfrm_cre
     *  number of fermion creation indices in the SQ operator product
     * @param nfrm_ann
     *  number of fermion annihilation indices in the SQ operator product
     * @param nbos_cre
     *  number of boson creation indices in the SQ operator product
     * @param nbos_ann
     *  number of boson annihilation indices in the SQ operator product
     * @return
     *  the excitation signature
     */
    static constexpr uint_t encode(uint_t nfrm_cre, uint_t nfrm_ann, uint_t nbos_cre, uint_t nbos_ann) {
        return (nfrm_cre > c_nop_mask_frm || nfrm_ann > c_nop_mask_frm ||
                nbos_cre > c_nop_mask_bos || nbos_ann > c_nop_mask_bos) ?
               ~0ul : nfrm_cre | (nfrm_ann << c_nbit_nop_frm) |
                      (nbos_cre << (2 * c_nbit_nop_frm)) |
                      (nbos_ann << (2 * c_nbit_nop_frm + c_nbit_nop_bos));
    }

public:
    constexpr OpSig(uint_t opsig): m_i(opsig){}

    constexpr OpSig(): OpSig(0ul){}

    constexpr OpSig(const OpSig& opsig): OpSig(opsig.m_i){}

    constexpr OpSig(OpSig&& opsig): OpSig(opsig.m_i){}

    constexpr OpSig(uintp_t frm, uintp_t bos): m_i(encode(frm.first, frm.second, bos.first, bos.second)){}

//    operator const uint_t& ()  const {
//        return m_i;
//    }
//    operator const bool ()  const {
//        return m_i < c_ndistinct;
//    }

    bool is_valid() const {
        return m_i < c_ndistinct;
    }

    /**
     * @return
     *  the number of fermion creation indices in the SQ operator product encoded within exsig
     */
    constexpr uint_t nfrm_cre() const {
        return c_nop_mask_frm & m_i;
    }

    /**
     * @return
     *  the number of fermion annihilation indices in the SQ operator product encoded within exsig
     */
    constexpr uint_t nfrm_ann() const {
        return c_nop_mask_frm & (m_i >> c_nbit_nop_frm);
    }

    /**
     * @return
     *  the number of boson creation indices in the SQ operator product encoded within exsig
     */
    constexpr uint_t nbos_cre() const {
        return c_nop_mask_bos & (m_i >> (2 * c_nbit_nop_frm));
    }

    /**
     * @return
     *  the number of boson annihilation indices in the SQ operator product encoded within exsig
     */
    constexpr uint_t nbos_ann() const {
        return c_nop_mask_bos & (m_i >> (2 * c_nbit_nop_frm + c_nbit_nop_bos));
    }

    /**
     * @return
     *  the total number of fermion indices in the SQ operator product encoded within exsig
     */
    constexpr uint_t nfrm() const {
        return nfrm_cre() + nfrm_ann();
    }

    /**
     * @return
     *  the total number of boson indices in the SQ operator product encoded within exsig
     */
    constexpr uint_t nbos() const {
        return nbos_cre() + nbos_ann();
    }

    /**
     * @return
     *  the total number of operators of any particle type in the SQ operator product encoded within exsig
     */
    constexpr uint_t nop() const {
        return nfrm() + nbos();
    }

    /**
     * @return
     *  true if the exsig has no boson operators
     */
    constexpr bool is_pure_frm() const {
        return !(nbos_cre() + nbos_ann());
    }

    /**
     * @return
     *  true if the exsig has no fermion operators
     */
    constexpr bool is_pure_bos() const {
        return !(nfrm_cre() + nfrm_ann());
    }

    /**
     * @return
     *  true if the exsig represents a fermion number-conserving operator product
     */
    constexpr bool conserves_nfrm() const {
        return nfrm_cre() == nfrm_ann();
    }

    /**
     * @return
     *  true if the exsig represents a boson number-conserving operator product
     */
    constexpr bool conserves_nbos() const {
        return nbos_cre() == nbos_ann();
    }

    /**
     * @return
     *  assuming this is a rank signature, the number of contributing promotions in the fermionic operators
     */
    constexpr uint_t ncontrib_frm() const {
        return integer::min(nfrm_cre(), nfrm_ann()) + 1;
    }

    /**
     * @return
     *  assuming this is a rank signature, the number of contributing promotions in the bosonic operators
     */
    constexpr uint_t ncontrib_bos() const {
        return integer::min(nbos_ann(), nbos_cre()) + 1;
    }

    /**
     * @return
     *  assuming this is a rank signature, the OpSig with the smallest nop which contributes (promotes to this ranksig)
     */
    OpSig base() const {
        return {
                {
                    nfrm_cre() - (ncontrib_frm() - 1),
                    nfrm_ann() - (ncontrib_frm() - 1)
                },
                {
                    nbos_cre() - (ncontrib_bos() - 1),
                    nbos_ann() - (ncontrib_bos() - 1)
                }
        };
    }

    /**
     * @param nfrm
     *  number of fermion operators to add (both cre and ann)
     * @param nbos
     *  number of boson operators to add (both cre and ann)
     * @return
     *  new OpSig with the additional operators
     */
    OpSig add_ops(uint_t nfrm, uint_t nbos) const {
        return {{nfrm_cre() + nfrm, nfrm_ann() + nfrm}, {nbos_cre() + nbos, nbos_ann() + nbos}};
    }

    constexpr bool contribs_to_frm(const OpSig& ranksig) const {
        return (nfrm_cre() <= ranksig.nfrm_cre()) && (
                (ranksig.nfrm_cre() - nfrm_cre()) == (ranksig.nfrm_ann() - nfrm_ann()));
    }

    constexpr bool contribs_to_bos(const OpSig& ranksig) const {
        return (nbos_cre() <= ranksig.nbos_cre()) && (
                (ranksig.nbos_cre() - nbos_cre()) == (ranksig.nbos_ann() - nbos_ann()));
    }

    constexpr bool contribs_to(const OpSig& ranksig) const {
        return contribs_to_frm(ranksig) && contribs_to_bos(ranksig);
    }

    constexpr bool takes_contribs_from(const OpSig& exsig) const {
        return exsig.contribs_to(*this);
    }

    /**
     * @return
     *  the OpSig representing the hermitian conjugate of the operator represented by this OpSig
     */
    OpSig conj() const {
        return {{nfrm_ann(), nfrm_cre()}, {nbos_cre(), nbos_ann()}};
    }

    str_t to_string() const {
        if (m_i >= c_ndistinct) return "invalid";
        return std::to_string(nfrm_cre()) + std::to_string(nfrm_ann()) +
               std::to_string(nbos_cre()) + std::to_string(nbos_ann());
    }
};

#endif //M7_OPSIG_H
