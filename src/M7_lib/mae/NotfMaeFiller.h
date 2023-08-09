//
// Created by anderson on 09/08/2023.
//

#ifndef M7_NOTFMAEFILLER_H
#define M7_NOTFMAEFILLER_H

#include "M7_lib/bilinear/Rdms.h"
#include "M7_lib/wavefunction/WalkerTable.h"
#include "M7_lib/util/BitsetIntersection.h"

class NotfMaeFiller {

    /**
     * Histogrammable set of determinants of which to compute the outer product in filling the MAEs
     */
    const Table<MbfWeightRow>& m_hist;
    /**
     * Wrapper around all RDM and MRPT2 intermediate objects filled by this class
     */
    Rdms* m_rdms;  // uninitialised pointer to RDM object?
    /**
     * offset from beginning of hist MBFs for the local MPI rank
     */
    const uint_t m_ind_displ;
    /**
     * number of hist MBFs for the local MPI rank
     */
    const uint_t m_ind_count;

    /**
     * bitsets encoding the indices of all MBFs which contain an electron in the spinorbital corresponding to the index
     * of the outer vector.
     * e.g. m_occ_bitsets[p] is a vector of uint_t ints which have set bits at the position indices where the MBF in
     * m_hist has spinorb p occupied
     */
    const v_t<uintv_t> m_occ_bitsets;
    /**
     * as above, but restricted to a rank-local region of the histogrammable set
     */
    const v_t<uintv_t> m_partial_occ_bitsets;

    /**
     * m_occ_bitsets expressed in isect form
     */
    const v_t<v_t<uintp_t>> m_occ_isects;  // these exclude the zero words, e.g. 0010 | 0011 | 0000 | 1010 -> [[0, 0010], [1, 0011], [3, 1010]]
    /**
     * m_partial_occ_bitsets expressed in isect form
     */
    const v_t<v_t<uintp_t>> m_partial_occ_isects;


    v_t<uintv_t> make_occ_bitsets(uint_t displ, uint_t count) const;  // declare existence of member function?

    v_t<uintv_t> make_occ_bitsets() const;

    /**
     * creates a hash table of the (nelec - rank)-electron determinants which resolve the identity between the creation
     * and annihilation SQ operators in the normal-ordered product
     * @param ann_ispinorbs
     *  strictly ascending order vector of spinorb indices in the annihilation operators taken only from MBFs in the
     *  index range [m_ind_displ, m_ind_displ + m_ind_count) of m_hist
     * @param ann_isect
     *  result of intersecting all bitsets corresponding to the ann_ispinorbs
     * @param cre_ispinorbs
     *  strictly ascending order vector of spinorb indices in the annihilation operators taken from the entirety of m_hist
     * @param cre_isect
     *  result of intersecting all bitsets corresponding to the cre_ispinorbs
     */
    void resolve_identity(const uintv_t& /*ann_ispinorbs*/, const v_t<uintp_t>& /*ann_isect*/,
                          const uintv_t& /*cre_ispinorbs*/, const v_t<uintp_t>& /*cre_isect*/) {
        REQUIRE_TRUE_ALL(m_rdms, "RDMs object must be non-null");
    }


    /**
     * enumerate all normal-ordered products of fermion spinorb SQ operators which contribute to the RDMs to be filled.
     * for each of these combinations, dispatch resolve identity which matches the creation and annihilation set intersections
      // what does the above line mean exactly? Does the definition of resolve_identity change depending on the input? It always receives the same types?
     * and fills the non-zero contributions to all RDMs
     */
    void fill() {
        // [&] => every argument to this lambda can be found by reference outside its scope
        auto fn = [&](const uintv_t& ao, const v_t<uintp_t>& ais, const uintv_t& co, const v_t<uintp_t>& cis) {
            resolve_identity(ao, ais, co, cis);
        };
        fill_foreach_isect_pair(fn);
    }

public:
    // what is this? Neither a variable nor a function nor an enum? A constructor?
    NotfMaeFiller(const Table<MbfWeightRow>& hist, Rdms* rdms=nullptr):
        m_hist(hist), m_rdms(rdms),
        m_ind_displ(mpi::evenly_shared_displ(m_hist.nrow_in_use())),
        m_ind_count(mpi::evenly_shared_count(m_hist.nrow_in_use())),
        m_occ_bitsets(make_occ_bitsets()),
        m_partial_occ_bitsets(make_occ_bitsets(m_ind_displ, m_ind_count)),
        m_occ_isects(bitset_isect::make_isects(m_occ_bitsets)),
        m_partial_occ_isects(bitset_isect::make_isects(m_partial_occ_bitsets)) {

    }

    // `static` avoids namespace conflict by restricting linkage to this translation unit
    static void fill(const Table<MbfWeightRow>& hist, Rdms* rdms=nullptr) {
        NotfMaeFiller filler(hist, rdms);  // initialiase a filler
        filler.fill();
    }

    template<typename fn_t>
    void fill_foreach_isect_pair(const fn_t& fn) {
        functor::assert_prototype<void(const uintv_t& /*ann_ispinorbs*/, const v_t<uintp_t>& /*ann_isect*/,
                                       const uintv_t& /*cre_ispinorbs*/, const v_t<uintp_t>& /*cre_isect*/)>(fn);

        // todo: Arta
        // uintv_t ann_ispinorbs;
        // v_t<uintp_t> ann_isect;
        // uintv_t cre_ispinorbs;
        // v_t<uintp_t> cre_isect;
        // some body of iterative or recursive code that enumerates and forms sets {
        //      fn(ann_ispinorbs, ann_isect, cre_ispinorbs, cre_isect);
        // }
    }

};


#endif //M7_NOTFMAEFILLER_H
