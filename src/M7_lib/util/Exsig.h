//
// Created by rja on 12/06/22.
//

#ifndef M7_UTIL_EXSIG_H
#define M7_UTIL_EXSIG_H

#include "utils.h"
#include "Integer.h"

namespace utils {

    /**
     * functions related to "excitation signatures" of connections between many-body basis functions
     */
    namespace exsig {

        using namespace defs;

        /**
         * compactly expresses an arbitrary SQ operator product as a single integer given some compile-time constant numbers
         * of bits for each element. e.g. if nbit_exsig_nop_frm = 3 and nbit_exsig_nop_bos = 1, then 2x3+2x1 = 8 bits are
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
        static constexpr size_t encode(size_t nfrm_cre, size_t nfrm_ann, size_t nbos_cre, size_t nbos_ann) {
            return (nfrm_cre > exsig_nop_mask_frm || nfrm_ann > exsig_nop_mask_frm ||
                    nbos_cre > exsig_nop_mask_bos || nbos_ann > exsig_nop_mask_bos) ?
                   ~0ul : nfrm_cre | (nfrm_ann << nbit_exsig_nop_frm) |
                          (nbos_cre << (2 * nbit_exsig_nop_frm)) |
                          (nbos_ann << (2 * nbit_exsig_nop_frm + nbit_exsig_nop_bos));
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  the number of fermion creation indices in the SQ operator product encoded within exsig
         */
        static constexpr size_t decode_nfrm_cre(size_t exsig) {
            return exsig_nop_mask_frm & exsig;
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  the number of fermion annihilation indices in the SQ operator product encoded within exsig
         */
        static constexpr size_t decode_nfrm_ann(size_t exsig) {
            return exsig_nop_mask_frm & (exsig >> nbit_exsig_nop_frm);
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  the number of boson creation indices in the SQ operator product encoded within exsig
         */
        static constexpr size_t decode_nbos_cre(size_t exsig) {
            return exsig_nop_mask_bos & (exsig >> (2 * nbit_exsig_nop_frm));
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  the number of boson annihilation indices in the SQ operator product encoded within exsig
         */
        static constexpr size_t decode_nbos_ann(size_t exsig) {
            return exsig_nop_mask_bos & (exsig >> (2 * nbit_exsig_nop_frm + nbit_exsig_nop_bos));
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  the total number of fermion indices in the SQ operator product encoded within exsig
         */
        static constexpr size_t decode_nfrm(size_t exsig) {
            return decode_nfrm_cre(exsig) + decode_nfrm_ann(exsig);
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  the total number of boson indices in the SQ operator product encoded within exsig
         */
        static constexpr size_t decode_nbos(size_t exsig) {
            return decode_nbos_cre(exsig) + decode_nbos_ann(exsig);
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  the total number of operators of any particle type in the SQ operator product encoded within exsig
         */
        static constexpr size_t decode_nop(size_t exsig) {
            return decode_nfrm(exsig) + decode_nbos(exsig);
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  true if the exsig has no boson operators
         */
        static constexpr bool is_pure_frm(size_t exsig) {
            return !(decode_nbos_cre(exsig) || decode_nbos_ann(exsig));
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  true if the exsig has no fermion operators
         */
        static constexpr bool is_pure_bos(size_t exsig) {
            return !(decode_nfrm_cre(exsig) || decode_nfrm_ann(exsig));
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  true if the exsig represents a fermion number-conserving operator product
         */
        static constexpr bool conserves_nfrm(size_t exsig) {
            return decode_nfrm_cre(exsig) == decode_nfrm_ann(exsig);
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  true if the exsig represents a boson number-conserving operator product
         */
        static constexpr bool conserves_nbos(size_t exsig) {
            return decode_nbos_cre(exsig) == decode_nbos_ann(exsig);
        }

        static constexpr size_t ncontrib_frm(size_t exsig) {
            return integer::min(decode_nfrm_cre(exsig), decode_nfrm_ann(exsig)) + 1;
        }

        static constexpr size_t ncontrib_bos(size_t exsig) {
            return integer::min(decode_nbos_cre(exsig), decode_nbos_ann(exsig)) + 1;
        }

        static constexpr size_t base_exsig(size_t exsig) {
            return encode(
                    decode_nfrm_cre(exsig) - (ncontrib_frm(exsig) - 1),
                    decode_nfrm_ann(exsig) - (ncontrib_frm(exsig) - 1),
                    decode_nbos_cre(exsig) - (ncontrib_bos(exsig) - 1),
                    decode_nbos_ann(exsig) - (ncontrib_bos(exsig) - 1)
            );
        }

        static constexpr size_t add_ops(size_t exsig, size_t nfrm, size_t nbos) {
            return encode(decode_nfrm_cre(exsig) + nfrm, decode_nfrm_ann(exsig) + nfrm,
                          decode_nbos_cre(exsig) + nbos, decode_nbos_ann(exsig) + nbos);
        }

        static constexpr bool contribs_to_frm(size_t exsig, size_t ranksig) {
            return (decode_nfrm_cre(exsig) <= decode_nfrm_cre(ranksig)) && (
                    decode_nfrm_cre(ranksig) - decode_nfrm_cre(exsig) == decode_nfrm_ann(ranksig) - decode_nfrm_ann(exsig));
        }

        static constexpr bool contribs_to_bos(size_t exsig, size_t ranksig) {
            return (decode_nbos_cre(exsig) <= decode_nbos_cre(ranksig)) && (
                    decode_nbos_cre(ranksig) - decode_nbos_cre(exsig) == decode_nbos_ann(ranksig) - decode_nbos_ann(exsig));
        }

        static constexpr bool contribs_to(size_t exsig, size_t ranksig) {
            return contribs_to_frm(exsig, ranksig) && contribs_to_bos(exsig, ranksig);
        }

        /**
         * @param exsig
         *  excitation signature
         * @return
         *  the exsig representing the hermitian conjugate of the operator represented by the argument
         */
        static constexpr size_t hermconj(size_t exsig) {
            return utils::exsig::encode(decode_nfrm_ann(exsig), decode_nfrm_cre(exsig), decode_nbos_ann(exsig),
                                       decode_nbos_cre(exsig));
        }

        static std::string to_string(size_t exsig) {
            if (exsig > nexsig) return "invalid";
            return std::to_string(decode_nfrm_cre(exsig)) + std::to_string(decode_nfrm_ann(exsig)) +
                   std::to_string(decode_nbos_cre(exsig)) + std::to_string(decode_nbos_ann(exsig));
        }

        static constexpr size_t ex_single = encode(1, 1, 0, 0);
        static constexpr size_t ex_double = encode(2, 2, 0, 0);
        static constexpr size_t ex_triple = encode(3, 3, 0, 0);
        static constexpr size_t ex_1100 = ex_single;
        static constexpr size_t ex_2200 = ex_double;
        static constexpr size_t ex_3300 = ex_triple;
        static constexpr size_t ex_1101 = encode(1, 1, 0, 1);
        static constexpr size_t ex_1110 = encode(1, 1, 1, 0);
        static constexpr size_t ex_0001 = encode(0, 0, 0, 1);
        static constexpr size_t ex_0010 = encode(0, 0, 1, 0);
        static constexpr size_t ex_0011 = encode(0, 0, 1, 1);
        static constexpr size_t ex_0022 = encode(0, 0, 2, 2);
    }
}

#endif //M7_UTIL_EXSIG_H