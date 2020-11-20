//
// Created by Robert John Anderson on 2020-04-11.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H

#include <src/core/sample/PRNG.h>
#include <src/core/excitgen/ExcitationGenerator.h>
#include <src/core/excitgen/UniformSingles.h>
#include <src/core/excitgen/HeatBathDoubles.h>
#include <src/core/sample/WeightedDrawer.h>
#include "Propagator.h"

class StochasticPropagator : public Propagator {

protected:
    PRNG m_prng;
    std::vector<std::unique_ptr<ExcitationGenerator>> m_exgens;
    const double &m_min_spawn_mag;
    std::unique_ptr<WeightedDrawer> m_exgen_drawer;

public:
    StochasticPropagator(const Hamiltonian &ham, const Options &opts) :
            Propagator(ham, opts), m_prng(opts.prng_seed, opts.prng_ngen),
            m_min_spawn_mag(opts.min_spawn_mag) {

        m_exgens.push_back(std::unique_ptr<ExcitationGenerator>(
                new UniformSingles(&m_ham, m_prng)));
        if (opts.excit_gen == "pchb") {
            m_exgens.push_back(std::unique_ptr<ExcitationGenerator>(
                    new HeatBathDoubles(&m_ham, m_prng)));
        }

        m_exgen_drawer = std::unique_ptr<WeightedDrawer>(new WeightedDrawer(m_exgens.size(), m_prng));

        const defs::prob_t prob_boson = 0.2; // TODO: make dynamic.
        if (m_exgens.size() == 2)
            m_exgen_drawer->set(m_magnitude_logger.m_psingle);
        else if (m_exgens.size() == 3)
            m_exgen_drawer->set(m_magnitude_logger.m_psingle, 1.0 - m_magnitude_logger.m_psingle - prob_boson);

        std::cout << "Excitation class probability breakdown " << utils::to_string(m_exgen_drawer->m_probs) << std::endl;
    }

    void diagonal(Wavefunction &m_wf, const size_t &irow) override {
        bool flag_deterministic = m_wf.m_walkers.m_flags.m_deterministic(irow);
        auto hdiag = m_wf.m_walkers.m_hdiag(irow);
        if (flag_deterministic) {
            m_wf.scale_weight(irow, 1 - (hdiag - m_shift) * tau());
        } else {
            // the probability that each unit walker will die
            auto death_rate = (hdiag - m_shift) * tau();
            if (death_rate < 0.0 || death_rate > 1.0) {
                // clone  / create antiparticles continuously
                m_wf.scale_weight(irow, 1 - death_rate);
            } else {
                auto weight = m_wf.m_walkers.m_weight(irow, 0, 0);
                // kill stochastically
                m_wf.set_weight(irow, m_prng.stochastic_round(weight, 1.0) * (1 - death_rate));
            }
        }
    }


    template<typename T>
    size_t get_nattempt(const T &weight) {
        static_assert(std::is_floating_point<T>::value, "template arg must be floating point");
        /*
         * We want to make nattempt = ceil(|weight|) spawning attempts.
         * can't rely on std::abs to provide the right answer in the case of complex arithmetic with
         * a real-valued FermionHamiltonian and an integral weight, since the sqrt function is not the exact
         * inverse of squaring in finite precision arithmetic!
         */
        if (weight < 0) return (-weight) < 1 ? 1 : std::round(-weight);
        else return (weight) < 1 ? 1 : std::round(weight);
    }

    template<typename T>
    size_t get_nattempt(const std::complex<T> &weight) {
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
        for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
            size_t iexgen = m_exgen_drawer->draw();
            auto &exgen = m_exgens[iexgen];
            if (!exgen->draw(src_onv, m_dst_onv, m_occ, m_vac, prob, helem, m_aconn)) continue;
            prob*=m_exgen_drawer->prob(iexgen);
            auto delta = -(weight / (defs::ham_comp_t) nattempt) * tau() * helem / prob;
            if (consts::float_is_zero(delta)) continue;
            m_wf.add_spawn(m_dst_onv, delta, flag_initiator, flag_deterministic);
        }
    }
};


//            switch (nexcit) {
//                case 1:
//                    valid = m_exgen->draw_single(src_onv, m_dst_onv, m_occ, m_vac, prob, helem, m_aconn);
//                    if (!valid) break;
//                    ASSERT(prob >= 0.0 && prob <= 1.0)
//                    prob *= m_magnitude_logger.m_psingle;
//                    ASSERT(!consts::float_nearly_zero(prob, 1e-14));
//                    break;
//                case 2:
//                    // TODO: don't need m_vac for doubles.
//                    valid = m_exgen->draw_double(src_onv, m_dst_onv, m_occ, prob, helem, m_aconn);
//                    if (!valid) break;
//                    ASSERT(prob >= 0.0 && prob <= 1.0)
//                    prob *= 1.0 - m_magnitude_logger.m_psingle;
//                    break;
//                default:
//                    throw std::runtime_error("invalid excitation rank");
//            }
//
//#ifdef VERBOSE_DEBUGGING
//            std::cout << consts::verb << consts::chevs << "EXCITATION GENERATED" << std::endl;
//        std::cout << consts::verb << "excitation rank:         " << nexcit << std::endl;
//        std::cout << consts::verb << "is valid:                " << string_utils::yn(valid) << std::endl;
//#endif
//
//            if (!valid) continue;
//            ASSERT(!consts::float_is_zero(prob))
//            auto delta = -(weight / (defs::ham_comp_t) nattempt) * tau() * helem / prob;
//#ifdef VERBOSE_DEBUGGING
//            std::cout << consts::verb << "probability:             " << prob << std::endl;
//        std::cout << consts::verb << "H matrix element:        " << helem << std::endl;
//        std::cout << consts::verb << "continuous delta:        " << delta << std::endl;
//#endif
//            delta = m_prng.stochastic_threshold(delta, m_min_spawn_mag);
//#ifdef VERBOSE_DEBUGGING
//            std::cout << consts::verb << "delta post-thresh:       " << delta << std::endl;
//#endif
//
//            ASSERT(consts::floats_equal(delta, -(weight / (defs::ham_comp_t) nattempt) * tau() * helem / prob)
//                   || consts::float_is_zero(delta) || consts::float_is_zero(delta - m_min_spawn_mag))
//
//            if (consts::float_is_zero(delta)) continue;
//            ASSERT(m_dst_onv.nsetbit() == src_onv.nsetbit())
//
//            m_wf.add_spawn(m_dst_onv, delta, flag_initiator, flag_deterministic);
//            m_magnitude_logger.log(nexcit, helem, prob);

#endif //M7_STOCHASTICPROPAGATOR_H