//
// Created by Robert J. Anderson on 03/06/2021.
//

#ifndef M7_SHIFT_H
#define M7_SHIFT_H

#include <queue>

#include <M7_lib/io/InteractiveVariable.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/wavefunction/Wavefunction.h>
#include <M7_lib/parallel/Epoch.h>

#include "MagnitudeLogger.h"
#include "M7_lib/wavefunction/Reference.h"
#include "M7_lib/wavefunction/Wavefunction.h"


namespace shift {
    /**
     * responsible for defining and updating the shift subtracted from H elements in the diagonal "death" step of propagation
     */
    struct ShiftSpace {
        /**
         * shift space index: 0 is the index of the most senior space - the one which determines when the variable mode
         * the epoch starts
         */
        const uint_t m_ispace;
        /**
         * update period (MC cycles)
         */
        const uint_t m_period;
        /**
         * the numbers of walkers on each WF part in the last period is stored so that the growth rate can be computed
         */
        buffered::Numbers<wf_comp_t, c_ndim_wf> m_nw_last_period;
        /**
         * values of the diagonal shift for each WF part
         */
        buffered::Numbers<ham_comp_t, c_ndim_wf> m_values;
        /**
         * every shift update protocol has an associated target number of walkers
         */
        const wf_comp_t m_nw_target;

        ShiftSpace(const NdFormat<c_ndim_wf>& wf_fmt, uint_t ispace, uint_t period, ham_comp_t init, wf_comp_t nw_target) :
            m_ispace(ispace), m_period(period), m_nw_last_period(wf_fmt.m_shape, std::numeric_limits<wf_comp_t>::max()),
            m_values(wf_fmt.m_shape, init), m_nw_target(nw_target) {
            m_nw_last_period.clear();
        }

        const ham_comp_t& operator[](uint_t ipart) const {
            return m_values[ipart];
        }

        bool is_period_cycle(uint_t icycle) const {
            return !(icycle % m_period);
        }

        virtual str_t enter_variable_mode_str(const wf::Vectors& wf, uint_t icycle, uint_t ipart) const = 0;

        /**
         * compute the change in all parts of the shift value
         * @param wf
         *  wavefunction whose population growth defines the change in shift
         * @param icycle
         *  MC cycle index
         * @param tau
         *  current timestep
         * @param variable_mode
         *  epochs begin when the shift value is to be modulated in order to satisfy some condition
         */
        virtual void update(const wf::Vectors& wf, uint_t icycle, double tau, Epochs& variable_mode) = 0;
    };

    struct GrowthBased : ShiftSpace {
        /**
         * damping factor in the shift update expression
         */
        const double m_damp_fac;
        /**
         * if using target-driven damping, this will be y^2/4 where y is m_damp_fac, else it will be 0
         */
        const double m_target_damp_fac;

        GrowthBased(const NdFormat<c_ndim_wf>& wf_fmt, uint_t ispace, uint_t period, ham_comp_t init,
                    wf_comp_t nw_target, double damp_fac, bool target_damp):
            ShiftSpace(wf_fmt, ispace, period, init, nw_target),
            m_damp_fac(damp_fac), m_target_damp_fac(target_damp ? math::pow<2>(m_damp_fac)/1.0 : 0.0){}

        str_t enter_variable_mode_str(const wf::Vectors& wf, uint_t icycle, uint_t ipart) const override {
            const auto nw = wf.m_stats.m_nw.total()[ipart];
            return logging::format("Variable shift triggered for WF part {}. Cycle {} nw: {}, cycle {} nw: {}",
                              ipart, icycle - 1, wf.m_stats.m_nw.prev_total()[ipart], icycle, nw);
        }

        void update(const wf::Vectors& wf, uint_t icycle, double tau, Epochs& variable_mode) override {
            for (uint_t ipart = 0ul; ipart < variable_mode.nelement(); ++ipart) {
                /*
                 * at the beginning of cycle i - where this update is performed, Nw_i is not available directly since the loop
                 * over the current occupied list has yet to be performed. Nw_i is required so compute S_i, so we must get it by
                 * adding the difference in Nw due to the application of cycle i-1 propagator.
                 */
                auto nw = wf.m_stats.m_nw.total()[ipart];
                /*
                 * number of cycles since last update
                 */
                uint_t a = 0ul;

                if (variable_mode[ipart].update(icycle, std::abs(nw) >= std::abs(m_nw_target))) {
                    if (icycle) {
                        logging::info("Variable shift triggered for WF part {}. Cycle {} nw: {}, cycle {} nw: {}",
                                      ipart, icycle - 1, wf.m_stats.m_nw.prev_total()[ipart], icycle, nw);
                    } else {
                        logging::info("Variable shift triggered immediately for WF part {}.", ipart);
                    }
                    a = icycle % m_period;
                }

                if (is_period_cycle(icycle)) a = m_period;

                if (variable_mode[ipart] && a) {
                    // if this is the first cycle, we have no growth rate, so set it to 1 i.e. "unchanged"
                    auto rate = icycle ? nw / m_nw_last_period[ipart] : 1.0;
                    m_values[ipart] -= m_damp_fac * std::log(std::abs(rate)) / (tau * a);
                    if (m_target_damp_fac != 0.0) {
                        rate = nw / m_nw_target;
                        m_values[ipart] -= m_target_damp_fac * std::log(std::abs(rate)) / (tau * a);
                    }
                }
            }
            if (is_period_cycle(icycle)) m_nw_last_period = wf.m_stats.m_nw.total();
        }
    };

    struct RefWeightFixing : ShiftSpace {

    };
}

/**
 * responsible for defining and updating the shift subtracted from the diagonal H elements in
 * propagation
 */
struct Shifts {
    shift::GrowthBased m_growth_based;
    /**
     * when this epoch begins, the shift is allowed to vary for the corresponding WF part
     */
    Epochs m_variable_mode;
    /**
     * values of the diagonal shift for each space and for each WF part
     */
    buffered::Numbers<ham_comp_t, 1+c_ndim_wf> m_values;
    /**
     * threshold for a shift_space > 0 MBF to be promoted to shift space 0
     */
    const ham_comp_t m_log_enhancement_promote_thresh;

    Shifts(const conf::Shift &opts, const NdFormat<c_ndim_wf>& wf_fmt):
        m_growth_based(wf_fmt, 0, opts.m_period, opts.m_init, opts.m_nw_targets.m_value[0], opts.m_damp, opts.m_target_damp),
        m_variable_mode("variable shift mode", wf_fmt.m_nelement, "WF part"),
        m_values(wf_fmt.add_major_dim(nspace(), "shift space")),
        m_log_enhancement_promote_thresh(opts.m_log_enhancement_promote_thresh){}

    const shift::ShiftSpace& operator[](const Walker& walker) const {
        (void) walker;
        return m_growth_based;
    }

    void update(const wf::Vectors& wf, uint_t icycle, double tau) {
        auto& first_shift_space = m_growth_based;
        // first update the variable mode epochs for all WF parts.
//        for (uint_t ipart=0ul; ipart < wf.m_format.m_nelement; ++ipart) {
//            if (m_variable_mode[ipart].update(icycle, first_shift_space.enter_variable_mode(wf, ipart))) {
//                if (!icycle) logging::info("Variable shift triggered immediately for WF part {}.", ipart);
//                else logging::info(first_shift_space.enter_variable_mode_str(wf, icycle, ipart));
//            }
//        }
        first_shift_space.update(wf, icycle, tau, m_variable_mode);
        auto& format = m_values.m_format;
        for (uint_t ipart=0ul; ipart < wf.m_format.m_nelement; ++ipart)
            m_values[format.combine<2>(0, ipart)] = first_shift_space.m_values[ipart];
    }

    uint_t nspace() const {
        return 1ul;//m_opts.m_nw_targets.m_value.size();
    }

    Shifts& operator=(const ham_comp_t& v) {
        m_growth_based.m_values = v;
        return *this;
    }

    Shifts& operator+=(const ham_comp_t& v) {
        m_growth_based.m_values += v;
        return *this;
    }

};

#endif //M7_SHIFT_H
