//
// Created by Robert J. Anderson on 10/08/2021.
//

#ifndef M7_FERMIONPROMOTER_H
#define M7_FERMIONPROMOTER_H

#include <M7_lib/connection/Connections.h>
#include <M7_lib/field/Fields.h>

/**
 * given a Connection, iterate over common operators
 *
 * Fermionic and bosonic promotion require sorting of indices, fermionic additionally requires the number of swaps
 * involved in the sort to be recorded in order to compute the phase of the promotion
 *
 * for a fermion connection of excitation level N contributing to a fermion RDM of rank R (> N), a contribution to the
 * RDM is due for each of the (N_site-N) choose (R-N) occupied indices in common between the contributing bra and ket.
 *
 * elements of the m_cre and m_ann members of both fields::FermionMevInds conn::FrmOnv alike must be stored in
 * ascending order. However, whereas the connection arrays must not have any elements in common, the MEV uintv_t in general
 * do have elements in common, this is brought about by promotion.
 *
 * e.g.
 *     01234 56789
 *    (01101,01011)   <- ket determinant
 *    (00111,11001)   <- bra determinant
 *    ( x   ,   x )
 *    ann: [1, 8]
 *    (   x ,x    )
 *    cre: [3, 5]
 *    (  x x, x  x)
 *    com: [2, 4, 6, 9]
 *    phase: -1
 *
 *    promotion to 3-body RDM element:
 *    enumerate common indices i and count number of set indices between i and the beginning of
 *    the string.
 *
 */
struct FermionPromoter {
    /**
     * number of common creation-annihilation operator pairs to insert into connection
     */
    const uint_t m_nop_insert;
    /**
     * total number of possible combinations
     */
    const uint_t m_ncomb;
    /**
     * enumeration of all possible combinations
     */
    std::vector<mev_ind_t> m_all_combs;

    FermionPromoter(uint_t ncom, uint_t nop_insert);

private:
    /**
     * @param icomb
     *  combination index
     * @return
     *  const pointer to beginning of combination in m_all_combs
     */
    const mev_ind_t *begin(const uint_t &icomb) const;

public:
    /**
     * apply the icomb-th promotion to the connection given and store the result uintv_t
     * @param icomb
     *  combination index
     * @param conn
     *  connection (no repeated SQ operator indices between ann and cre vectors)
     * @param com
     *  indices in common between bra and ket fermion ONVs
     * @param frm_inds
     *  MEV index field
     * @return
     *  antisymmetric phase associated with sorting both ann and cre to ascending order
     */
    bool apply(const uint_t &icomb, const conn::FrmOnv &conn,
               const FrmOps& com, MaeIndsPair &frm_inds) const;

};


#endif //M7_FERMIONPROMOTER_H
