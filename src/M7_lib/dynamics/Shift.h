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

        ShiftSpace(const NdFormat<c_ndim_wf>& wf_fmt, uint_t ispace, uint_t period, ham_comp_t init, wf_comp_t nw_target);
        virtual ~ShiftSpace() = default;

        const ham_comp_t& operator[](uint_t ipart) const {
            return m_values[ipart];
        }

        bool is_period_cycle(uint_t icycle) const {
            return !(icycle % m_period);
        }

        uint_t ncycle_this_update(uint_t ipart, uint_t icycle, const Epochs& variable_mode) const;

        virtual void update_variable_mode(const wf::Vectors& wf, uint_t icycle, Epochs& variable_mode) = 0;

        void update(const wf::Vectors& wf, uint_t icycle, double tau, const Epochs& variable_mode);
    protected:
        /**
         * compute the change in one part of the shift value
         * @param wf
         *  wavefunction whose population growth defines the change in shift
         * @param ipart
         *  WF part index
         * @param icycle
         *  MC cycle index
         * @param tau
         *  current timestep
         * @param variable_mode
         *  epochs begin when the shift value is to be modulated in order to satisfy some condition
         */
        virtual void update_part(const wf::Vectors& wf, uint_t ipart, uint_t icycle, double tau, const Epochs& variable_mode) = 0;

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
                    wf_comp_t nw_target, double damp_fac, bool target_damp);

        void update_variable_mode(const wf::Vectors& wf, uint_t icycle, Epochs& variable_mode) override;

    protected:
        void update_part(const wf::Vectors& wf, uint_t ipart, uint_t icycle, double tau, const Epochs& variable_mode) override;
    };

    struct RefWeightFixing : ShiftSpace {
        RefWeightFixing(const NdFormat<c_ndim_wf>& wf_fmt, uint_t ispace, uint_t period, ham_comp_t init, wf_comp_t nw_target);

        void update_variable_mode(const wf::Vectors& wf, uint_t icycle, Epochs& variable_mode) override;

    protected:
        void update_part(const wf::Vectors& wf, uint_t ipart, uint_t icycle, double, const Epochs& variable_mode) override;
    };

    struct ValueFixing : ShiftSpace {
        ValueFixing(const NdFormat<c_ndim_wf>& wf_fmt, uint_t ispace, uint_t period, ham_comp_t init);

        void update_variable_mode(const wf::Vectors& wf, uint_t icycle, Epochs& variable_mode) override;

    protected:
        void update_part(const wf::Vectors&, uint_t, uint_t, double, const Epochs&) override {}
    };
}

/**
 * responsible for defining and updating the shift subtracted from the diagonal H elements in propagation
 */
struct Shifts {
    v_t<std::unique_ptr<shift::ShiftSpace>> m_spaces;
    /**
     * when this epoch begins, the shift is allowed to vary for the corresponding WF part
     */
    Epochs m_variable_mode;
    /**
     * values of the diagonal shift for each space and for each WF part
     */
    buffered::Numbers<ham_comp_t, 1+c_ndim_wf> m_values;
    /**
     * if enhancement takes the form exp(a*delta S), then a is the quantity in the range [0, 1] that damps it.
     * if a = 0, no enhancement is done
     * id a = 1, full enhancement is done
     */
    const ham_comp_t m_enhancement_damp;
    /**
     * if a log enhancement factor exceeds this value, the walker is promoted to the s0 space
     */
    const ham_comp_t m_s0_promote_thresh;
    /**
     * minimum shift values
     */
    const v_t<ham_comp_t> m_floors;

    Shifts(const conf::Shift &opts, const NdFormat<c_ndim_wf>& wf_fmt);

    const shift::ShiftSpace* operator[](const Walker& walker) const;

    void update(const wf::Vectors& wf, uint_t icycle, double tau);

    Shifts& operator=(const ham_comp_t& v);

    Shifts& operator+=(const ham_comp_t& v);

    void promote_to_s0_if_high_enhancement(wf::Vectors& wf, Walker& walker) const;

};

#endif //M7_SHIFT_H
