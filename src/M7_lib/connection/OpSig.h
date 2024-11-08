//
// Created by rja on 26/10/22.
//

#ifndef M7_OPSIG_H
#define M7_OPSIG_H

#include "M7_lib/util/Integer.h"


#if 0
struct OpSign {
    uint_t m_i;
    uint_t m_nfrm;
    uint_t m_npos;
    uint_t m_nbos_cre;
    uint_t m_nbos_ann;

};


namespace opsign {
    static constexpr uint_t c_nsig = 15;

    static constexpr uint_t c_1f = 1ul;
    static constexpr uint_t c_2f = 2ul;
    static constexpr uint_t c_3f = 3ul;
    static constexpr uint_t c_4f = 4ul;
    static constexpr uint_t c_11b = 5ul;
    static constexpr uint_t c_22b = 6ul;
    static constexpr uint_t c_01b = 7ul;
    static constexpr uint_t c_10b = 8ul;
    static constexpr uint_t c_1f_01b = 9ul;
    static constexpr uint_t c_1f_10b = 10ul;
    static constexpr uint_t c_1p = 11ul;
    static constexpr uint_t c_1f_1p = 12ul;
    static constexpr uint_t c_2f_1p = 13ul;

    static constexpr uint_t c_sing = c_1f;
    static constexpr uint_t c_doub = c_2f;
    static constexpr uint_t c_trip = c_3f;
    static constexpr uint_t c_quad = c_4f;

    //                                                        0  1  2  3  4  5  6  7  8  9  10 11 12 13 14
    static constexpr std::array<uint_t, c_nsig> c_nfrm     = {0, 1, 2, 3, 4, 0, 0, 0, 0, 1, 1, 0, 1, 2, ~0ul};
    static constexpr std::array<uint_t, c_nsig> c_npos     = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, ~0ul};
    static constexpr std::array<uint_t, c_nsig> c_nbos_cre = {0, 0, 0, 0, 0, 1, 2, 0, 1, 0, 1, 0, 0, 0, ~0ul};
    static constexpr std::array<uint_t, c_nsig> c_nbos_ann = {0, 0, 0, 0, 0, 1, 2, 1, 0, 1, 0, 0, 0, 0, ~0ul};

    static constexpr std::array<bool, c_nsig>   c_pure_frm = {1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    static constexpr std::array<bool, c_nsig>   c_pure_pos = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
    static constexpr std::array<bool, c_nsig>   c_pure_bos = {1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
    static constexpr std::array<uint_t, c_nsig> c_iconj    = {0, 1, 2, 3, 4, 0, 0, 0, 0, 1, 1, 0, 1, 2, ~0ul};

    static constexpr uint_t c_ncontribs_to_ranksigs = 38;
    static constexpr std::array<uint_t, c_ncontribs_to_ranksigs> c_ranksigs = { 
        // exsig                                          end offset
        /*  0 */   0,  1,  2,  3,  4,  5,  6, 11, 12, 13,  // 10
        /*  1 */   1,  2,  3,  4, 12, 13,                  // 16
        /*  2 */   2,  3,  4, 13,                          // 20
        /*  3 */   3,  4,                                  // 22
        /*  4 */   4,                                      // 23
        /*  5 */   5,  6,                                  // 25
        /*  6 */   6,                                      // 26
        /*  7 */   7,  9,                                  // 28
        /*  8 */   8,  10,                                 // 30
        /*  9 */   9,                                      // 31
        /* 10 */  10,                                      // 32
        /* 11 */  11, 12, 13,                              // 35
        /* 12 */  12, 13,                                  // 37
        /* 13 */  13                                       // 38
        /* 14 */                                           // 38
    };
    static constexpr std::array<uint_t, c_nsig> c_ranksig_end_offsets = {10, 16, 20, 22, 23, 25, 26, 28, 30, 31, 32, 35, 37, 38, 38};

    const uint_t* ranksig_cbegin(uint_t exsig) {
        return c_ranksigs.data() + (exsig ? c_ranksig_end_offsets[exsig-1] : 0ul);
    }

    const uint_t* ranksig_cend(uint_t exsig) {
        return c_ranksigs.data() + c_ranksig_end_offsets[exsig];
    }

    constexpr uint_t nranksig(uint_t exsig) {
        return c_ranksig_end_offsets[exsig] - (exsig ? c_ranksig_end_offsets[exsig-1] : 0ul);
    }

}
#endif

namespace opsig {
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
}

/**
 * OpSigs encode the rank or excitation level of a second-quantized operator string in a compact form. The exact number
 * and range of possible single integer representations of OpSigs is determined by the numbers of bits defined in the
 * above namespace. The representable ranks are therefore limited to those useful in ascertaining connectivity with
 * respect to a particular operator e.g. Hamiltonian or RDM
 */
class OpSig {
    uint_t m_i;

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
        return (nfrm_cre > opsig::c_nop_mask_frm || nfrm_ann > opsig::c_nop_mask_frm ||
                nbos_cre > opsig::c_nop_mask_bos || nbos_ann > opsig::c_nop_mask_bos) ?
               ~0ul : nfrm_cre | (nfrm_ann << opsig::c_nbit_nop_frm) |
                      (nbos_cre << (2 * opsig::c_nbit_nop_frm)) |
                      (nbos_ann << (2 * opsig::c_nbit_nop_frm + opsig::c_nbit_nop_bos));
    }

public:

    struct Pair {
        uint_t m_ncre, m_nann;
    };
    constexpr explicit OpSig(uint_t opsig): m_i(opsig){}

    constexpr OpSig(): OpSig(0ul){}

    constexpr OpSig(Pair frm, Pair bos): m_i(encode(frm.m_ncre, frm.m_nann, bos.m_ncre, bos.m_nann)){}

    constexpr OpSig(const OpSig& opsig): OpSig(opsig.m_i){}

    OpSig& operator=(const OpSig& other) {
        return (m_i = other.m_i, *this);
    }

    constexpr OpSig(OpSig&& opsig): OpSig(opsig.m_i){}

    OpSig& operator=(OpSig&& other) noexcept {
        return (m_i = other.m_i, *this);
    }

    bool operator==(const OpSig& other) const {
        return other.m_i == m_i;
    }

    bool operator!=(const OpSig& other) const {
        return other.m_i != m_i;
    }

    constexpr operator const uint_t& () const {
        return m_i;
    }

    /**
     * @return
     *  the number of fermion creation indices in the SQ operator product encoded within exsig
     */
    constexpr uint_t nfrm_cre() const {
        return opsig::c_nop_mask_frm & m_i;
    }

    /**
     * @return
     *  the number of fermion annihilation indices in the SQ operator product encoded within exsig
     */
    constexpr uint_t nfrm_ann() const {
        return opsig::c_nop_mask_frm & (m_i >> opsig::c_nbit_nop_frm);
    }

    /**
     * @return
     *  the number of boson creation indices in the SQ operator product encoded within exsig
     */
    constexpr uint_t nbos_cre() const {
        return opsig::c_nop_mask_bos & (m_i >> (2 * opsig::c_nbit_nop_frm));
    }

    /**
     * @return
     *  the number of boson annihilation indices in the SQ operator product encoded within exsig
     */
    constexpr uint_t nbos_ann() const {
        return opsig::c_nop_mask_bos & (m_i >> (2 * opsig::c_nbit_nop_frm + opsig::c_nbit_nop_bos));
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
    constexpr OpSig base() const {
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
    constexpr OpSig add_ops(uint_t nfrm, uint_t nbos) const {
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
    constexpr OpSig conj() const {
        return {{nfrm_ann(), nfrm_cre()}, {nbos_ann(), nbos_cre()}};
    }

    str_t to_string() const {
        if (m_i >= opsig::c_ndistinct) return "invalid";
        return std::to_string(nfrm_cre()) + std::to_string(nfrm_ann()) +
               std::to_string(nbos_cre()) + std::to_string(nbos_ann());
    }

    template<typename fn_t>
    static void foreach(const fn_t& fn) {
        functor::assert_prototype<void(const OpSig&)>(fn);
        for (uint_t ifrm_cre = 0ul; ifrm_cre <= opsig::c_nop_mask_frm; ++ifrm_cre) {
            for (uint_t ifrm_ann = 0ul; ifrm_ann <= opsig::c_nop_mask_frm; ++ifrm_ann) {
                for (uint_t ibos_cre = 0ul; ibos_cre <= opsig::c_nop_mask_bos; ++ibos_cre) {
                    for (uint_t ibos_ann = 0ul; ibos_ann <= opsig::c_nop_mask_bos; ++ibos_ann) {
                        fn(OpSig({ifrm_cre, ifrm_ann}, {ibos_cre, ibos_ann}));
                    }
                }
            }
        }
    }
};

namespace opsig {
    static constexpr OpSig c_0000 ({0, 0}, {0, 0});
    static constexpr OpSig c_invalid (~0ul);

    static constexpr OpSig c_1100 ({1, 1}, {0, 0});
    static constexpr OpSig c_2200 ({2, 2}, {0, 0});
    static constexpr OpSig c_3300 ({3, 3}, {0, 0});
    static constexpr OpSig c_4400 ({4, 4}, {0, 0});

    static constexpr OpSig c_zero = c_0000;
    static constexpr OpSig c_sing = c_1100;
    static constexpr OpSig c_doub = c_2200;
    static constexpr OpSig c_trip = c_3300;
    static constexpr OpSig c_quad = c_4400;

    static constexpr OpSig frm(uint_t ncre, uint_t nann) {
        return {{ncre, nann}, {0ul, 0ul}};
    }
    static constexpr OpSig frm(uint_t nop) {
        return frm(nop, nop);
    }
    static constexpr OpSig bos(uint_t ncre, uint_t nann) {
        return {{0ul, 0ul}, {ncre, nann}};
    }
    static constexpr OpSig bos(uint_t nop) {
        return bos(nop, nop);
    }

    static constexpr OpSig c_1101 ({1, 1}, {0, 1});
    static constexpr OpSig c_1110 ({1, 1}, {1, 0});
    static constexpr OpSig c_0001 ({0, 0}, {0, 1});
    static constexpr OpSig c_0010 ({0, 0}, {1, 0});
    static constexpr OpSig c_0011 ({0, 0}, {1, 1});
    static constexpr OpSig c_0022 ({0, 0}, {2, 2});
}

#endif //M7_OPSIG_H
