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
    Rdms* m_rdms = nullptr;
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
     * m_occ_bitsets expressed in siv form
     */
    const v_t<v_t<uintp_t>> m_occ_sivs;
    /**
     * m_partial_occ_bitsets expressed in siv form
     */
    const v_t<v_t<uintp_t>> m_partial_occ_sivs;
    /**
     * working object for MBF connections
     */
    conn::Mbf m_work_conn;

    /**
     * hash table used in the resolution of the identity in the (nelectron - rank) electron Hilbert space
     */
    buffered::MappedTable<MbfWeightRow> m_ri_map;

    v_t<uintv_t> make_occ_bitsets(uint_t displ, uint_t count) const;

    v_t<uintv_t> make_occ_bitsets() const;

    void refresh_ann_map(const uintv_t& ann_ispinorbs, const v_t<uintp_t>& ann_siv, field::RdmInds& rdm_inds);

    void probe_ann_map(const uintv_t& cre_ispinorbs, const v_t<uintp_t>& cre_siv, field::RdmInds& rdm_inds);

    /**
     * creates a hash table of the (nelec - rank)-electron determinants due to the annihilation operators acting on the
     * ket hist vector, then probes this table with the result of the creation operators acting on the bra which
     * resolves the identity between the creation and annihilation SQ operators in the normal-ordered product
     * @param ann_ispinorbs
     *  strictly ascending order vector of spinorb indices in the annihilation operators taken only from MBFs in the
     *  index range [m_ind_displ, m_ind_displ + m_ind_count) of m_hist
     * @param ann_siv
     *  result of intersecting all bitsets corresponding to the ann_ispinorbs
     * @param cre_ispinorbs
     *  strictly ascending order vector of spinorb indices in the annihilation operators taken from the entirety of m_hist
     * @param cre_siv
     *  result of intersecting all bitsets corresponding to the cre_ispinorbs
     * @param rdm_inds
     *  result of intersecting all bitsets corresponding to the cre_ispinorbs
     */
    void resolve_identity(const uintv_t& ann_ispinorbs, const v_t<uintp_t>& ann_siv,
                          const uintv_t& cre_ispinorbs, const v_t<uintp_t>& cre_siv, field::RdmInds& rdm_inds);

    /**
     * enumerate all normal-ordered products of fermion spinorb SQ operators which contribute to the RDMs to be filled.
     * for each of these combinations, dispatch resolve identity which matches the creation and annihilation set intersections
     * and fills the non-zero contributions to all RDMs
     */
    void fill();

    wf_comp_t get_norm() const;

public:
    explicit NotfMaeFiller(const Table<MbfWeightRow>& hist, Rdms* rdms=nullptr);

    static void fill(const Table<MbfWeightRow>& hist, Rdms* rdms=nullptr);

    template<typename fn_t>
    void fill_foreach_set_pair(const fn_t& fn, const v_t<OpSig>& rdm_opsigs) {
        functor::assert_prototype<void(const uintv_t &ann_ispinorbs, const v_t<uintp_t> &ann_siv,
                                       const uintv_t &cre_ispinorbs, const v_t<uintp_t> &cre_siv,
                                       field::RdmInds& rdm_inds)>(fn);
        using namespace bitset_isect;
        for (const auto& opsig: rdm_opsigs) {
            if (!opsig) continue;
            buffered::RdmInds rdm_inds(opsig);
            const auto rank = opsig.nfrm_cre();
            auto ann_fn = [&](const uintv_t &ann_ispinorbs, const siv_t &ann_siv) -> void {
                auto cre_fn = [&](const uintv_t &cre_ispinorbs, const siv_t &cre_siv) -> void {
                    fn(ann_ispinorbs, ann_siv, cre_ispinorbs, cre_siv, rdm_inds);
                };
                foreach_unique(m_occ_bitsets, m_occ_sivs, rank, true, cre_fn);
            };
            foreach_unique(m_occ_bitsets, m_partial_occ_sivs, rank, true, ann_fn);
        }
    }

    /**
     * compute the Fermi phase associated with a "half excitation" i.e. either the creation or annihilation strings.
     * Regardless of whether the operators are creation or annihilation in the RDM matrix element definition, they
     * are always applied to histogrammed determinants as annihilation operators
     * @param ispinorbs
     *  indices of the spin orbitals in the given MBF
     * @param ihist_mbf
     *  integer index of the MBF as a row in the m_hist table.
     * @return
     *  true if the Fermi phase of the "half excitation" is -1
     */
    bool half_excit_phase(const uintv_t& ispinorbs, const field::Mbf& mbf);
    /**
     * @param ihist_mbf
     *  integer index of the MBF as a row in the m_hist table.
     */
    bool half_excit_phase(const uintv_t& ispinorbs, uint_t ihist_mbf);
};


#endif //M7_NOTFMAEFILLER_H
