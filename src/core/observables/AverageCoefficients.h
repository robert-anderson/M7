//
// Created by rja on 03/03/2021.
//

#ifndef M7_AVERAGECOEFFICIENTS_H
#define M7_AVERAGECOEFFICIENTS_H

#include <src/core/basis/Connections.h>
#include <src/core/config/FciqmcConfig.h>
#include "src/core/field/Fields.h"
#include "src/core/table/BufferedTable.h"
#include "src/core/table/BufferedFields.h"
#include "MevTable.h"

struct AvCoeffsOneExlvl : BufferedTable<MevRow<defs::wf_t>, true> {
    buffered::FermionMevInds m_working_inds;
    defs::wf_t m_ref_coeff = 0.0;

    AvCoeffsOneExlvl(size_t nann, size_t ncre, size_t nvalue, size_t nbucket = 100) :
            BufferedTable<MevRow<defs::wf_t>, true>("average coefficients", {{nann, ncre, nvalue}, nbucket}),
            m_working_inds({nann, ncre}) {
        REQUIRE_EQ_ALL(nann, ncre, "different creation and annihilation operator numbers not currently supported");
    }

    LookupResult operator[](const conn::Basic<0> &key) {
        set_working_inds(key);
        return MappedTable<MevRow<defs::wf_t>>::operator[](m_working_inds);
    }

    size_t insert(const conn::Basic<0> &key) {
        set_working_inds(key);
        return MappedTable<MevRow<defs::wf_t>>::insert(m_working_inds);
    }

    void add(const conn::Basic<0> &key, const defs::wf_t& contrib) {
        auto irow = *(*this)[key];
        if (irow==~0ul) irow = insert(key);
        m_row.jump(irow);
        m_row.m_values[0]+=contrib;
    }

    size_t nop() const {
        return m_working_inds.m_ann.m_size;
    }

    std::vector<std::string> h5_field_names() const {
        return {m_row.m_inds.m_ann.m_name,
                m_row.m_inds.m_cre.m_name,
                m_row.m_values.m_name};
    }

    void save(hdf5::GroupWriter& gw) const {
        this->write(gw, std::to_string(nop()), h5_field_names());
    }

private:
    void set_working_inds(const conn::Basic<0> &key) {
        m_working_inds = {key.ann(), key.cre()};
    }
};


struct AverageCoefficients {
    const fciqmc_config::Observables& m_opts;
    const size_t m_max_exlvl;
    std::vector<AvCoeffsOneExlvl> m_av_coeffs;
    AverageCoefficients(const fciqmc_config::Observables& opts) :
        m_opts(opts), m_max_exlvl(opts.m_av_coeffs.m_max_exlvl){
        m_av_coeffs.reserve(m_max_exlvl+1);
        for (size_t i=0ul; i<m_max_exlvl; ++i) m_av_coeffs.emplace_back(m_max_exlvl, m_max_exlvl, 1);
    }

    void add(const conn::Basic<0> &key, const defs::wf_t& contrib) {
        auto nop = key.nexcit();
        m_av_coeffs[nop].add(key, contrib);
    }

    void save(std::string fname="av_coeffs.h5") const {
        hdf5::FileWriter fw(fname);
        hdf5::GroupWriter gw("exlvls", fw);
        for (auto& it: m_av_coeffs) it.save(gw);
    }

    void save(size_t i) const {
        save(log::format("av_coeffs.{}.h5", i));
    }

    bool all_stores_empty() const {
        for (const auto& it: m_av_coeffs) if(it.m_hwm) return false;
        return true;
    };
};

#endif //M7_AVERAGECOEFFICIENTS_H
