//
// Created by Robert J. Anderson on 03/06/2021.
//

#include "Shift.h"
#include "M7_lib/util/Math.h"

shift::ShiftSpace::ShiftSpace(const NdFormat<c_ndim_wf>& wf_fmt, uint_t ispace, uint_t period, ham_comp_t init,
                              wf_comp_t nw_target) :
        m_ispace(ispace), m_period(period), m_nw_last_period(wf_fmt.m_shape, std::numeric_limits<wf_comp_t>::max()),
        m_values(wf_fmt.m_shape, init), m_nw_target(nw_target) {
    m_nw_last_period.clear();
}

shift::GrowthBased::GrowthBased(const NdFormat<c_ndim_wf>& wf_fmt, uint_t ispace, uint_t period, ham_comp_t init,
        wf_comp_t nw_target, double damp_fac, bool target_damp) : ShiftSpace(wf_fmt, ispace, period, init, nw_target),
        m_damp_fac(damp_fac), m_target_damp_fac(target_damp ? math::pow<2>(m_damp_fac)/1.0 : 0.0){}


uint_t shift::ShiftSpace::ncycle_this_update(uint_t ipart, uint_t icycle, const Epochs& variable_mode) const {
    if (!variable_mode[ipart]) return 0;
    /*
     * number of cycles since last update
     */
    uint_t a = 0ul;
    /*
     * perform update for full period if icycle is an integer multiple of the period, else check if there is an
     * initial update for a partial period
     */
    if (is_period_cycle(icycle)) a = m_period;
    else if (variable_mode[ipart].started_this_cycle(icycle)) a = icycle % m_period;
    return a;
}

uint_t get_iflat(const wf::Vectors& wf, uint_t ipart, uint_t ispace) {
    const auto& format = wf.m_stats.m_nw_by_shift_space.m_format;
    return format.combine<2>(ipart, ispace);
}

void shift::ShiftSpace::update(const wf::Vectors& wf, uint_t icycle, double tau, const Epochs& variable_mode) {
    for (uint_t ipart=0ul; ipart < variable_mode.nelement(); ++ipart){
        update_part(wf, ipart, icycle, tau, variable_mode);
        if (is_period_cycle(icycle))
            m_nw_last_period[ipart] = wf.m_stats.m_nw_by_shift_space.total()[get_iflat(wf, ipart, m_ispace)];
    }
}

void shift::GrowthBased::update_variable_mode(const wf::Vectors& wf, uint_t icycle, Epochs& variable_mode) {
    for (uint_t ipart = 0ul; ipart < variable_mode.nelement(); ++ipart) {
        const auto nw = wf.m_stats.m_nw_by_shift_space.total()[get_iflat(wf, ipart, m_ispace)];
        if (variable_mode[ipart].update(icycle, std::abs(nw) >= std::abs(m_nw_target))) {
            if (icycle) {
                logging::info("Variable shift triggered for WF part {}. Cycle {} nw: {}, cycle {} nw: {}",
                              ipart, icycle - 1, wf.m_stats.m_nw.prev_total()[ipart], icycle, nw);
            } else {
                logging::info("Variable shift triggered immediately for WF part {}.", ipart);
            }
        }
    }
}

void shift::GrowthBased::update_part(const wf::Vectors& wf, uint_t ipart, uint_t icycle, double tau, const Epochs& variable_mode) {
    const auto nw = wf.m_stats.m_nw_by_shift_space.total()[get_iflat(wf, ipart, m_ispace)];
    const auto a = ncycle_this_update(ipart, icycle, variable_mode);
    if (a) {
        // if this is the first cycle, or we otherwise have inf growth rate, so set it to 1 i.e. "unchanged"
        auto rate = std::abs(m_nw_last_period[ipart]) == 0.0 ? 1.0 : nw / m_nw_last_period[ipart];
        m_values[ipart] -= m_damp_fac * std::log(std::abs(rate)) / (tau * a);
        if (m_target_damp_fac != 0.0) {
            rate = nw / m_nw_target;
            m_values[ipart] -= m_target_damp_fac * std::log(std::abs(rate)) / (tau * a);
        }
    }
}

shift::RefWeightFixing::RefWeightFixing(const NdFormat<c_ndim_wf>& wf_fmt, uint_t ispace, uint_t period,
                                        ham_comp_t init, wf_comp_t nw_target) :
        ShiftSpace(wf_fmt, ispace, period, init, nw_target) {
    REQUIRE_EQ_ALL(period, 1ul, "reference weight fixing shift protocol requires update period of 1 cycle");
}

void shift::RefWeightFixing::update_variable_mode(const wf::Vectors& wf, uint_t icycle, Epochs& variable_mode) {
    for (uint_t ipart = 0ul; ipart < variable_mode.nelement(); ++ipart) {
        const auto mag = std::abs(wf.m_refs[ipart].weight());
        if (variable_mode[ipart].update(icycle, mag >= m_nw_target)) {
            if (icycle) {
                logging::info("Variable shift triggered for WF part {}. Cycle {} ref magnitude: {}",
                              ipart, icycle, mag);
            }
            else logging::info("Variable shift triggered immediately for WF part {}.", ipart);
        }
    }
}

void shift::RefWeightFixing::update_part(const wf::Vectors& wf, uint_t ipart, uint_t icycle, double,
                                         const Epochs& variable_mode) {
    const auto a = ncycle_this_update(ipart, icycle, variable_mode);
    if (a) {
        const auto e = wf.reference_projected_energy(ipart);
        m_values[ipart] = e;
    }
}

shift::ValueFixing::ValueFixing(const NdFormat<c_ndim_wf>& wf_fmt, uint_t ispace, uint_t period, ham_comp_t init) :
        ShiftSpace(wf_fmt, ispace, period, init, 0.0){}

void shift::ValueFixing::update_variable_mode(const wf::Vectors&, uint_t icycle, Epochs& variable_mode) {
    for (uint_t ipart = 0ul; ipart < variable_mode.nelement(); ++ipart) {
        if (variable_mode[ipart].update(icycle, true))
            logging::info("Fixing shift immediately in shift space {} for WF part {}.", m_ispace, ipart);
    }
}

Shifts::Shifts(const conf::Shift& opts, const NdFormat<c_ndim_wf>& wf_fmt) :
        m_variable_mode("variable shift mode", wf_fmt.m_nelement, "WF part"),
        m_values(wf_fmt.add_major_dim(opts.m_nw_targets.m_value.size(), "shift space")),
        m_enhancement_damp(opts.m_enhancement_damp), m_s0_promote_thresh(opts.m_s0_promote_thresh){
    // if the first space is of the "fix ref weight" or "fix s0" type, add it explicitly
    if (opts.m_fix_s0)
        m_spaces.emplace_back(new shift::ValueFixing(wf_fmt, 0, opts.m_period, opts.m_init));
    else if (opts.m_fix_ref_weight)
        m_spaces.emplace_back(new shift::RefWeightFixing(wf_fmt, 0, opts.m_period, opts.m_init, opts.m_nw_targets.m_value[0]));

    uint_t ispace = m_spaces.size();

    // add all remaining spaces as growth-based shifts
    for (; ispace<opts.m_nw_targets.m_value.size(); ++ispace)
        m_spaces.emplace_back(
            new shift::GrowthBased(wf_fmt, ispace, opts.m_period,
            opts.m_init, opts.m_nw_targets.m_value[ispace], opts.m_damp, opts.m_target_damp));

    logging::info("Initialized {}", string::plural("shift space", m_spaces.size()));
}

const shift::ShiftSpace* Shifts::operator[](const Walker& walker) const {
    return m_spaces[walker.m_shift_space].get();
}

void Shifts::update(const wf::Vectors& wf, uint_t icycle, double tau) {
    // only S0 determines when variable shift mode begins
    m_spaces[0]->update_variable_mode(wf, icycle, m_variable_mode);
    for (uint_t ispace=0ul; ispace < m_spaces.size(); ++ispace) {
        m_spaces[ispace]->update(wf, icycle, tau, m_variable_mode);
        auto& format = m_values.m_format;
        for (uint_t ipart=0ul; ipart < wf.m_format.m_nelement; ++ipart) {
            DEBUG_ASSERT_FALSE(math::is_nan_or_inf(std::abs(m_spaces[ispace]->m_values[ipart])), "new shift is invalid");
            // constrain shift values relative to ispace 0
            m_spaces[ispace]->m_values[ipart] = std::min(m_spaces[ispace]->m_values[ipart], m_spaces[0]->m_values[ipart]);
            // constrain shift values above S1 to be exactly S0
            if (ispace > 1) m_spaces[ispace]->m_values[ipart] = m_spaces[0]->m_values[ipart];
            m_values[format.combine<2>(ispace, ipart)] = m_spaces[ispace]->m_values[ipart];
        }
    }
}

Shifts& Shifts::operator=(const ham_comp_t& v) {
    for (auto& ptr: m_spaces) ptr->m_values = v;
    return *this;
}

Shifts& Shifts::operator+=(const ham_comp_t& v) {
    for (auto& ptr: m_spaces) ptr->m_values += v;
    return *this;
}

void Shifts::promote_to_s0_if_high_enhancement(wf::Vectors& wf, Walker& walker) const {
    // zero thresh signifies no promotion allowed
    if (m_s0_promote_thresh == 0.0) return;
    if (walker.m_shift_space > 0 && walker.m_log_enhancement_fac >= m_s0_promote_thresh) {
        const auto fac = std::exp(m_enhancement_damp * walker.m_log_enhancement_fac);
        wf.scale_weight(walker, 0, fac, 0);
        walker.m_log_enhancement_fac = 0.0;
    }
}