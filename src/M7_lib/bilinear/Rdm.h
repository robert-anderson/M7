//
// Created by Robert J. Anderson on 11/08/2021.
//

#ifndef M7_RDM_H
#define M7_RDM_H

#include <M7_lib/basis/Suites.h>
#include <M7_lib/wavefunction/Reference.h>
#include <M7_lib/propagator/Propagator.h>
#include <M7_lib/mae/MaeTable.h>
#include <M7_lib/field/Fields.h>
#include <M7_lib/io/Archivable.h>
#include <M7_lib/hdf5/Node.h>

#include "FermionPromoter.h"

class Rdm : public communicator::MappedSend<MaeRow, MaeRow> {
public:
    /**
     * RDM accumulation requires knowledge of the size of the basis and the numbers of particles
     */
    const sys::Sector m_sector;
    /**
     * rank signature of the RDM, along with convenient decoded numbers of operators in the ranksig
     */
    const OpSig m_ranksig;
    const uint_t m_nfrm_cre, m_nfrm_ann, m_nbos_cre, m_nbos_ann;
    /**
     * index signature (different to ranksig if RDM is stored as an on-the-fly contraction)
     */
    const OpSig m_indsig;
    const uint_t m_nfrm_cre_ind, m_nfrm_ann_ind, m_nbos_cre_ind, m_nbos_ann_ind;

protected:
    bool m_ordered_inds = true;
    /**
     * enumerators of the promotions of contributing excitation signatures to the ranksig of the RDM
     */
    v_t<FermionPromoter> m_frm_promoters;
    /**
     * indices of the full second quantised string
     */
    buffered::MaeInds m_full_inds;
    /**
     * indices of any contracted intermediate, not the full second quantised string
     */
    buffered::MaeInds m_uncontracted_inds;
    const str_t m_name;

    /**
     * objects for table traversal
     */
    MaeRow m_send_row;
    MaeRow m_recv_row;
    MaeRow m_store_row;

    static uint_t nrow_estimate(uint_t nfrm_cre, uint_t nfrm_ann, uint_t nbos_cre, uint_t nbos_ann, sys::Size basis_size);

    static uint_t nrow_estimate(OpSig exsig, sys::Size basis_size);

    str_t name(str_t str, OpSig ranksig) const;

    void add_to_send_table(const field::MaeInds& inds, wf_t contrib);

    virtual void frm_make_contribs(const field::FrmOnv& /*src_onv*/, const conn::FrmOnv& /*conn*/,
                                   const com_ops::Frm& /*com*/, wf_t /*contrib*/) {}

    virtual void frmbos_make_contribs(const field::FrmBosOnv& /*src_onv*/, const conn::FrmBosOnv& /*conn*/,
                                      const com_ops::FrmBos& /*com*/, wf_t /*contrib*/) {}

    virtual void bos_make_contribs(const field::BosOnv& /*src_onv*/, const conn::BosOnv& /*conn*/,
                                   const com_ops::Bos& /*com*/, wf_t /*contrib*/) {}

public:

    str_t name() const {
        return name(m_name, m_ranksig);
    }

    Rdm(OpSig ranksig, OpSig indsig, sys::Sector sector, uint_t nvalue,
        DistribOptions dist_opts, Sizing store_sizing, Sizing comm_sizing, str_t name="");

    /**
     * @param opts
     *  RDM section of the config document
     * @param ranksig
     *  rank of the SQ operators in each contribution
     * @param sector
     *  dimensions of the stored basis and number of particles to use in enforcing probability-conserving trace
     * @param nvalue
     *  number of values to encode in each RDM element
     * @param name
     *  string identifier for logging and archive output. if empty, this is generated from the ranksig
     * @param indsig
     *  number of each species of SQ operator to store in the structure (equal to ranksig for ordinary, uncontracted RDMs)
     */
    Rdm(const conf::Rdms& opts, OpSig ranksig, OpSig indsig, sys::Sector sector, uint_t nvalue, str_t name="");

    void end_cycle();

    void save(hdf5::NodeWriter& gw) const;

    void make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn, const com_ops::Frm& com, wf_t contrib) {
        frm_make_contribs(src_onv, conn, com, contrib);
    }

    void make_contribs(const field::FrmBosOnv& src_onv, const conn::FrmBosOnv& conn,
                       const com_ops::FrmBos& com, wf_t contrib) {
        frmbos_make_contribs(src_onv, conn, com, contrib);
    }

    void make_contribs(const field::BosOnv& src_onv, const conn::BosOnv& conn,
                       const com_ops::Bos& com, wf_t contrib) {
        bos_make_contribs(src_onv, conn, com, contrib);
    }
};

class PureRdm : public Rdm {
public:
    /**
     * @param opts
     *  RDM section of the config document
     * @param ranksig
     *  rank of the SQ operators in each contribution
     * @param sector
     *  dimensions of the stored basis and number of particles to use in enforcing probability-conserving trace
     * @param nvalue
     *  number of values to encode in each RDM element
     * @param name
     *  string identifier for logging and archive output. if empty, this is generated from the ranksig
     * @param indsig
     *  number of each species of SQ operator to store in the structure (equal to ranksig for ordinary, uncontracted RDMs)
     */
    PureRdm(const conf::Rdms& opts, OpSig ranksig, sys::Sector sector, uint_t nvalue, str_t name=""):
            Rdm(opts, ranksig, ranksig, sector, nvalue, name){}

protected:
    void frm_make_contribs(const FrmOnv& src_onv, const conn::FrmOnv& conn,
                           const com_ops::Frm& com, wf_t contrib) override;

    void frmbos_make_contribs(const FrmBosOnv& src_onv, const conn::FrmBosOnv& conn,
                              const com_ops::FrmBos& com, wf_t contrib) override;

    void bos_make_contribs(const BosOnv& /*onv*/, const conn::BosOnv& /*bosOnv*/,
                           const com_ops::Bos& /*bos*/, wf_t /*wf*/) override {
        ABORT("not yet implemented");
    }
};

/**
 * To reduce the memory storage requirement of high-rank RDMs, they may be contracted with a tensor on the fly
 *  https://aip.scitation.org/doi/10.1063/1.5140086
 * such that instead of storing G[a, b, c,..., i, j, k,...], we instead store for example:
 * FG[b, c, ..., j, k, ...] = sum_ai F[a, i] G[a, b, c,..., i, j, k,...]
 * where F is a coefficient array of rank 1100
 */
class ContractedRdm : public Rdm {
public:

    /**
     * the maximum exsig contributing to this ContractedRdm
     */
    const OpSig m_max_contrib_exsig;

    /**
     * @param opts
     *  RDM section of the config document
     * @param ranksig
     *  rank of the SQ operators in each contribution
     * @param indsig
     *  rank of the SQ operators retained for element indexing after contraction
     * @param max_contrib_exsig
     *  highest exsig which can make a non-zero contribution
     * @param sector
     *  dimensions of the stored basis and number of particles to use in enforcing probability-conserving trace
     * @param nvalue
     *  number of values to encode in each RDM element
     * @param name
     *  string identifier for logging and archive output. if empty, this is generated from the ranksig
     * @param indsig
     *  number of each species of SQ operator to store in the structure (equal to ranksig for ordinary, uncontracted RDMs)
     */
    ContractedRdm(const conf::Rdms& opts, OpSig ranksig, OpSig indsig, OpSig max_contrib_exsig,
                  sys::Sector sector, uint_t nvalue, str_t name=""):
            Rdm(opts, ranksig, indsig, sector, nvalue, name), m_max_contrib_exsig(max_contrib_exsig){
        REQUIRE_TRUE(m_max_contrib_exsig.contribs_to(m_ranksig), "max contributing exsig is incompatible with ranksig");
    }

protected:

    void frmbos_make_contribs(const FrmBosOnv& src_onv, const conn::FrmBosOnv& conn,
                              const com_ops::FrmBos& com, wf_t contrib) override {
        // so far, there is no use for contracted mixed fermion-boson RDMs, delegate to fermion sector method
        make_contribs(src_onv.m_frm, conn.m_frm, com.m_frm, contrib);
    }

    void bos_make_contribs(const BosOnv& /*onv*/, const conn::BosOnv& /*bosOnv*/,
                           const com_ops::Bos& /*bos*/, wf_t /*wf*/) override {
        // so far, there is no use for contracted boson RDMs
    }
};

#endif //M7_RDM_H
