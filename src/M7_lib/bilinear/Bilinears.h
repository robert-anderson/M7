//
// Created by Robert J. Anderson on 11/08/2021.
//

#ifndef M7_BILINEARS_H
#define M7_BILINEARS_H

#include <utility>

#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/field/Fields.h>

#include "Rdms.h"
#include "SpectralMoment.h"

/**
 * Projection onto a trial wavefunction is sufficient for the estimation of many-body expectation values if the operator
 * in question commutes with the Hamiltonian. When the operator does not commute with the Hamiltonian, the projection
 * by a trial wavefunction is just an approximation called a "mixed estimator".
 *
 * In this general case, the only way to obtain estimates which are exact in the infinite sampling limit is to use the
 * wavefunction itself as the trial wavefunction, and thus we obtain estimates of a quantity which is "bilinear" in the
 * wavefunction.
 *
 * This is statistically problematic when the WF is stochastically propagated, since products of correlated walker
 * populations introduce a systematic bias. This is overcome by replicating walker populations for each root.
 *
 * Two different types of multidimensional bilinear estimator are currently implemented.
 *  1. reduced density matrices (RDMs)
 *  2. spectral moments (SpecMoms)
 *
 * In M7, RDMs are defined as the expectation value of a normal ordered product of second quantized (SQ) operators and
 * may include fermionic spin orbital indices or bosonic modes as free indices.
 *
 * SpecMoms are defined as matrix elements of some power of the Hamiltonian between perturbed wavefunctions, where the
 * perturbing operators have arbitrary rank but are limited to fermionic spin orbital indices.
 */
struct Bilinears {
    Rdms m_rdms;
    SpecMoms m_spec_moms;
    buffered::Numbers<wf_t, c_ndim_root> m_total_norm;
    //const bool m_explicit_hf_conns;

    //std::array<std::unique_ptr<SpectralMoment>, c_ndistinct> m_specmoms;
    //const uintv_t m_specmom_exsigs;

private:
    /**
     * @param string
     *  excitation signature or fermion operator rank defined in the configuration document
     * @return
     *  integer excitation signature
     */
    static uint_t parse_exsig(const str_t &string) {
        REQUIRE_TRUE_ALL(string.size() == 1 || string.size() == 4, "invalid exsig string specification");
        if (string.size() == 1) {
            uint_t rank = string::parse_decimal_digit(string.c_str());
            REQUIRE_LE_ALL(rank, exsig::c_nop_mask_frm, "number of fermion operators exceeds limit");
            return exsig::encode(rank, rank, 0, 0);
        }
        uint_t nfrm_cre = string::parse_decimal_digit(string.c_str());
        REQUIRE_LE_ALL(nfrm_cre, exsig::c_nop_mask_frm, "number of fermion creation operators exceeds limit");
        uint_t nfrm_ann = string::parse_decimal_digit(string.c_str()+1);
        REQUIRE_LE_ALL(nfrm_ann, exsig::c_nop_mask_frm, "number of fermion annihilation operators exceeds limit");
        uint_t nbos_cre = string::parse_decimal_digit(string.c_str()+2);
        REQUIRE_LE_ALL(nbos_cre, exsig::c_nop_mask_bos, "number of boson creation operators exceeds limit");
        uint_t nbos_ann = string::parse_decimal_digit(string.c_str()+3);
        REQUIRE_LE_ALL(nbos_ann, exsig::c_nop_mask_bos, "number of boson annihilation operators exceeds limit");
        return exsig::encode(nfrm_cre, nfrm_ann, nbos_cre, nbos_ann);
    }
    /**
     * @param strings
     *  excitation signatures or fermion operator ranks defined in the configuration document
     * @return
     *  integer excitation signatures
     */
    static uintv_t parse_exsigs(const strv_t &strings) {
        uintv_t out;
        for (auto &string: strings) out.push_back(parse_exsig(string));
        DEBUG_ASSERT_EQ(out.size(), strings.size(),
                        "output should have the same number of exsigs as specification");
        return out;
    }

public:

    operator bool() const {
        return m_rdms || m_spec_moms;
    }

    /**
     * replication is needed if we are estimating any bilinears stochastically
     * @param stoch
     *  true if the wavefunction is propagated by a stochastic process
     * @return
     *  true if replication of the walker populations is a statistical necessity
     */
    bool need_replication(bool stoch) const {
        if (!stoch) return false;
        return m_rdms || m_spec_moms;
    }

    /**
     * parents only need to be sent in WF spawning table if we're going to be accumulating RDMs at some point
     */
    bool need_send_parents() const {
        return m_rdms;
    }

    Bilinears(const conf::Mae &opts, uintv_t rdm_ranksigs, uintv_t /*specmom_exsigs*/,
              sys::Sector sector, const Epoch &epoch) :
            m_rdms(opts.m_rdm, rdm_ranksigs, sector, epoch), m_spec_moms(opts.m_spec_mom), m_total_norm({1}) {}

    Bilinears(const conf::Mae &opts, sys::Sector sector, const Epoch &epoch) :
            Bilinears(opts, parse_exsigs(opts.m_rdm.m_ranks),
                      parse_exsigs(opts.m_spec_mom.m_ranks), sector, epoch) {}

    bool all_stores_empty() const {
        return m_rdms.all_stores_empty() && m_spec_moms.all_stores_empty();
    }

    void end_cycle() {
        m_rdms.end_cycle();
    }

    /**
     * Make all "diagonal" contributions
     * should only be called if:
     *  the WF row is about to be removed
     *  the end of a sampling period has been reached
     *  the end of the calculation has been reached (all WF rows are about to be removed)
     * @param mbf
     */
    void make_contribs(const field::Mbf &onv, const wf_t &ci_product) {
        m_total_norm[0] += ci_product;
        m_rdms.make_contribs(onv, onv, ci_product);
    }

    ham_comp_t estimate_energy(const Hamiltonian &ham) {
        return m_rdms.get_energy(ham);
    }
};

#endif //M7_BILINEARS_H
