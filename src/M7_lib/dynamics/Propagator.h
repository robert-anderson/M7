//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_PROPAGATOR_H
#define M7_PROPAGATOR_H

#include <M7_lib/conf/Conf.h>
#include <M7_lib/io/Archivable.h>
#include <M7_lib/basis/Suites.h>

#include "Shift.h"

/**
 * This class is not optionally archivable,
 * if an archive is being saved or checkpointed, the propagator data must be included
 * if an archive is being loaded, the propagator data must be read and verified against the data on file to ensure
 * compatibility with any objects reinstated from disk
 */
class Propagator : public Archivable {
protected:
    double m_tau;
public:
    const NdFormat<defs::ndim_wf> m_wf_fmt;
    const Hamiltonian &m_ham;
    Shift m_shift;
    const sys::Sector m_sector;
    /**
     * exponent in the Gutzwiller-like importance sampling
     */
    const double m_imp_samp_exp;
    /*
     * working objects
     */
    mutable suite::Mbfs m_dst;
    mutable suite::Conns m_conn;

    Propagator(const conf::Document &opts, const Hamiltonian &ham, const Wavefunction &wf) :
            Archivable("propagator", opts.m_archive),
            m_tau(opts.m_propagator.m_tau_init),
            m_wf_fmt(wf.m_format),
            m_ham(ham),
            m_shift(opts, wf.m_format),
            m_sector(wf.m_sector),
            m_imp_samp_exp(opts.m_propagator.m_imp_samp_exp),
            m_dst(m_sector),
            m_conn(ham.m_basis.size()){}

    virtual ~Propagator() {}

    virtual void diagonal(Wavefunction &wf, const size_t &ipart) = 0;

    virtual void off_diagonal(Wavefunction &wf, const size_t &ipart) = 0;

    virtual defs::ham_t round(const defs::ham_t &weight) {
        return weight;
    }

    const double &tau() const {
        return m_tau;
    }

    virtual void update(const size_t &icycle, const Wavefunction &wf);

    virtual size_t ncase_excit_gen() const {
        return 0;
    }

    virtual std::vector<defs::prob_t> excit_gen_case_probs() const {
        return {};
    }

    /**
     * do Gutzwiller-like importance sampling (if enabled) of the spawning probability "delta"
     * @param delta
     *  reference to the normal spawning probability
     * @param dst_mbf
     *  destination MBF
     * @param src_energy
     *  energy of the source MBF is in practice cached in the wavefunction table, so no need to recompute
     */
    void imp_samp_delta(defs::wf_t& delta, const field::Mbf& dst_mbf, defs::ham_comp_t src_energy) const;

    void imp_samp_delta(defs::wf_t& delta, const field::Mbf& src_mbf, const field::Mbf& dst_mbf) const;

private:
    void load_fn(hdf5::GroupReader &parent) override;

    void save_fn(hdf5::GroupWriter &parent) override;

};

#endif //M7_PROPAGATOR_H
