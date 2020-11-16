//
// Created by Robert John Anderson on 2020-04-11.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H

#include <src/core/sample/PRNG.h>
#include <src/core/sample/ExcitationGenerator.h>
#include "Propagator.h"

class StochasticPropagator : public Propagator {
    PRNG m_prng;
    std::unique_ptr<ExcitationGenerator> m_exgen = nullptr;
    const double &m_min_spawn_mag;

public:
    StochasticPropagator(const Hamiltonian &ham, const Options& opts):
    Propagator(ham, opts), m_prng(opts.prng_seed, opts.prng_ngen), m_min_spawn_mag(opts.min_spawn_mag){}

    void diagonal(Wavefunction &m_wf, const size_t &irow) override{
        bool flag_deterministic = m_wf.m_walkers.m_flags.m_deterministic(irow);
        auto hdiag = m_wf.m_walkers.m_hdiag(irow);
        if (flag_deterministic) {
            m_wf.scale_weight(irow, 1 - (hdiag - m_shift) * tau());
        }
        else {
            // the probability that each unit walker will die
            auto death_rate = (hdiag - m_shift) * tau();
            if (death_rate < 0.0 || death_rate > 1.0) {
                // clone  / create antiparticles continuously
                m_wf.scale_weight(irow, 1 - death_rate);
            }
            else {
                auto weight = m_wf.m_walkers.m_weight(irow, 0, 0);
                // kill stochastically
                m_wf.set_weight(irow, m_prng.stochastic_round(weight, 1.0) * (1 - death_rate));
            }
        }
    }


    template<typename T>
    size_t get_nattempt(const T& weight){
        static_assert(std::is_floating_point<T>::value, "template arg must be floating point");
        /*
         * We want to make nattempt = ceil(|weight|) spawning attempts.
         * can't rely on std::abs to provide the right answer in the case of complex arithmetic with
         * a real-valued FermionHamiltonian and an integral weight, since the sqrt function is not the exact
         * inverse of squaring in finite precision arithmetic!
         */
        if (weight<0) return (-weight) < 1 ? 1 : std::round(-weight);
        else return (weight) < 1 ? 1 : std::round(weight);
    }

    template<typename T>
    size_t get_nattempt(const std::complex<T>& weight) {
        if (m_ham.complex_valued()) return get_nattempt(std::abs(weight));
        else return get_nattempt(consts::real(weight));
    }

    void off_diagonal(Wavefunction &m_wf, const size_t &irow) override {
        auto weight = m_wf.m_walkers.m_weight(irow, 0, 0);
        ASSERT(!consts::float_is_zero(weight));
        ASSERT(consts::imag(weight) == 0.0 || m_ham.complex_valued())
        auto src_onv = m_wf.m_walkers.m_onv(irow);
        bool flag_initiator = m_wf.m_walkers.m_flags.m_initiator(irow, 0, 0);
        bool flag_deterministic = m_wf.m_walkers.m_flags.m_deterministic(irow);

        m_occ.update(src_onv);
        m_vac.update(src_onv);
        size_t nattempt = get_nattempt(weight);
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << "spawn attempts: " << nattempt << std::endl;
#endif
        defs::prob_t prob;
        defs::ham_t helem;
        bool valid = false;
        for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
            size_t nexcit = 2 - m_prng.stochastic_round(m_magnitude_logger.m_psingle, 1);
            switch (nexcit) {
                case 1:
                    valid = m_exgen->draw_single(src_onv, m_dst_onv, m_occ, m_vac, prob, helem, m_aconn);
                    if (!valid) break;
                    ASSERT(prob >= 0.0 && prob <= 1.0)
                    prob *= m_magnitude_logger.m_psingle;
                    ASSERT(!consts::float_nearly_zero(prob, 1e-14));
                    break;
                case 2:
                    // TODO: don't need m_vac for doubles.
                    valid = m_exgen->draw_double(src_onv, m_dst_onv, m_occ, prob, helem, m_aconn);
                    if (!valid) break;
                    ASSERT(prob >= 0.0 && prob <= 1.0)
                    prob *= 1.0 - m_magnitude_logger.m_psingle;
                    break;
                default:
                    throw std::runtime_error("invalid excitation rank");
            }

#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << consts::chevs << "EXCITATION GENERATED" << std::endl;
        std::cout << consts::verb << "excitation rank:         " << nexcit << std::endl;
        std::cout << consts::verb << "is valid:                " << string_utils::yn(valid) << std::endl;
#endif

            if (!valid) continue;
            ASSERT(!consts::float_is_zero(prob))
            auto delta = -(weight / (defs::ham_comp_t) nattempt) * tau() * helem / prob;
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << "probability:             " << prob << std::endl;
        std::cout << consts::verb << "H matrix element:        " << helem << std::endl;
        std::cout << consts::verb << "continuous delta:        " << delta << std::endl;
#endif
            delta = m_prng.stochastic_threshold(delta, m_min_spawn_mag);
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << "delta post-thresh:       " << delta << std::endl;
#endif

            ASSERT(consts::floats_equal(delta, -(weight / (defs::ham_comp_t) nattempt) * tau() * helem / prob)
                   || consts::float_is_zero(delta) || consts::float_is_zero(delta - m_min_spawn_mag))

            if (consts::float_is_zero(delta)) continue;
            ASSERT(m_dst_onv.nsetbit() == src_onv.nsetbit())

            m_wf.add_spawn(m_dst_onv, delta, flag_initiator, flag_deterministic);
            m_magnitude_logger.log(nexcit, helem, prob);
        }
    }

};


#if 0

#include <src/core/sample/PRNG.h>
#include <src/core/sample/UniformSampler.h>
#include <src/core/pchb/HeatBathSamplers.h>
#include "Propagator.h"
#include "FciqmcCalculation.h"


class StochasticPropagator : public Propagator {
    const double &m_min_spawn_mag;
public:
    PRNG m_prng;
    std::unique_ptr<ExcitationGenerator> m_exgen = nullptr;

    StochasticPropagator(FciqmcCalculation *fciqmc) :
            Propagator(fciqmc), m_min_spawn_mag(m_input.min_spawn_mag),
            m_prng(m_input.prng_seed, m_input.prng_ngen) {
        std::cout << "Initializing stochastic propagator" << std::endl;
        if (m_input.excit_gen=="pchb"){
            m_exgen = std::unique_ptr<ExcitationGenerator>(new HeatBathSamplers(m_ham.get(), m_prng));
        }
        else if (m_input.excit_gen=="uniform"){
            //m_exgen = std::unique_ptr<ExcitationGenerator>(new UniformSampler(m_ham.get(), m_prng));
        }
        else {
            throw std::runtime_error("invalid excit_gen specified");
        }

    }


    void off_diagonal(const DeterminantElement &src_det, const NumericElement<defs::ham_t> &weight,
                      SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) override;

    void diagonal(const NumericElement<defs::ham_comp_t> &hdiag, NumericElement<defs::ham_t> &weight,
                  bool flag_deterministic,
                  defs::ham_comp_t &delta_square_norm, defs::ham_comp_t &delta_nw) override;
};

#endif //M7_STOCHASTICPROPAGATOR_H
#endif //M7_STOCHASTICPROPAGATOR_H
