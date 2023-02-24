//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_PROPAGATOR_H
#define M7_PROPAGATOR_H

#include <M7_lib/conf/Conf.h>
#include <M7_lib/io/Archivable.h>
#include <M7_lib/basis/Suites.h>

#include "M7_lib/dynamics/Shift.h"
#include "Guide.h"
#include "M7_lib/wavefunction/Reference.h"
#include "M7_lib/wavefunction/Wavefunction.h"

/**
 * This class is not optionally archivable,
 * if an archive is being saved or checkpointed, the propagator data must be included
 * if an archive is being loaded, the propagator data must be read and verified against the data on file to ensure
 * compatibility with any objects reinstated from disk
 */
class Propagator {
protected:
    double m_tau;
public:
    const wf_comp_t m_nadd_initiator;
    const NdFormat<c_ndim_wf> m_wf_fmt;
    const Hamiltonian &m_ham;
    Shift m_shift;
    const sys::Sector m_sector;
    /**
     * guiding wavefunction used in the importance sampling
     */
    std::unique_ptr<guide::Wavefunction> m_imp_samp_guide;
    /*
     * working objects
     */
    mutable suite::Mbfs m_dst;
    mutable suite::Conns m_conn;

    std::unique_ptr<guide::Wavefunction> make_imp_samp_guide(const conf::GuideWavefunction& opts) const;

    Propagator(const conf::Document &opts, const Hamiltonian &ham, const wf::Fci &wf) :
            m_tau(opts.m_propagator.m_tau_init),
            m_nadd_initiator(opts.m_propagator.m_nadd),
            m_wf_fmt(wf.m_format),
            m_ham(ham),
            m_shift(opts, wf.m_format),
            m_sector(wf.m_sector),
            m_imp_samp_guide(make_imp_samp_guide(opts.m_propagator.m_imp_samp_guide)),
            m_dst(m_sector),
            m_conn(ham.m_basis.size()) {}

    virtual ~Propagator() {}

    virtual void diagonal(wf::Fci &wf, Walker& walker, const uint_t &ipart) = 0;

    virtual void off_diagonal(wf::Fci &wf, const Walker& walker, const uint_t &ipart) = 0;

    virtual ham_t round(const ham_t &weight) {
        return weight;
    }

    const double &tau() const {
        return m_tau;
    }

    virtual void update(uint_t icycle, const wf::Fci &wf, const wf::Refs& refs);

    virtual uint_t ncase_excit_gen() const {
        return 0;
    }

    virtual v_t<prob_t> excit_gen_case_probs() const {
        return {};
    }

    /**
     * use a guiding wavefunction to perform importance sampling
     * @param delta
     *  generated spawn contribution
     * @param src_ovlp
     *  pre-computed overlap between the guiding WF and the source basis function
     * @param dst_mbf
     *  generated destination basis function
     */
    void imp_samp_delta(wf_t &delta, ham_t src_ovlp, const field::Mbf &dst_mbf) const;

    void imp_samp_delta(wf_t &delta, const field::Mbf &src_mbf, const field::Mbf &dst_mbf) const;

    virtual hash::digest_t checksum_() const {return 0;}
    hash::digest_t checksum() const;

};

#endif //M7_PROPAGATOR_H
