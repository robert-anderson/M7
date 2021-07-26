//
// Created by rja on 03/03/2021.
//

#ifndef M7_REFEXCITS_H
#define M7_REFEXCITS_H

#include <src/core/connection/Connections.h>
#include <src/core/config/FciqmcConfig.h>
#include <src/core/io/Archivable.h>
#include "src/core/field/Fields.h"
#include "src/core/table/BufferedTable.h"
#include "src/core/table/BufferedFields.h"
#include "MevTable.h"

struct RefExcitsOneExlvl : BufferedTable<MevRow<defs::wf_t>, true> {
    /**
     * work space for converting between indexing vector type and the stored key type of the MEV tables
     */
    buffered::FermionMevInds m_working_inds;

    RefExcitsOneExlvl(size_t nann, size_t ncre, size_t nvalue, size_t nbucket = 100) :
            BufferedTable<MevRow<defs::wf_t>, true>("average coefficients", {{nann, ncre, nvalue}, nbucket}),
            m_working_inds({nann, ncre}) {
        REQUIRE_EQ_ALL(nann, ncre, "different creation and annihilation operator numbers not currently supported");
    }

    LookupResult operator[](const conn::FrmOnv &key) {
        set_working_inds(key);
        return MappedTable<MevRow<defs::wf_t>>::operator[](m_working_inds);
    }

    size_t insert(const conn::FrmOnv &key) {
        set_working_inds(key);
        return MappedTable<MevRow<defs::wf_t>>::insert(m_working_inds);
    }

    size_t nop() const {
        return m_working_inds.m_ann.m_size;
    }

    std::vector<std::string> h5_field_names() const {
        return {m_row.m_inds.m_ann.m_name,
                m_row.m_inds.m_cre.m_name,
                m_row.m_values.m_name};
    }

    using Table<MevRow<defs::wf_t>>::save;
    void save(hdf5::GroupWriter& gw) const {
        Table<MevRow<defs::wf_t>>::save(gw, std::to_string(nop()), h5_field_names());
    }

    void make_contribs(const conn::FrmOnv& conn, const defs::wf_t& contrib, const size_t& ipart) {
        auto irow = *(*this)[conn];
        if (irow==~0ul) irow = insert(conn);
        m_row.jump(irow);
        m_row.m_values[ipart]+=contrib;
    }

private:
    void set_working_inds(const conn::FrmOnv &key) {
        m_working_inds = {key.ann(), key.cre()};
    }
};


struct RefExcits : Archivable {
    const fciqmc_config::RefExcits& m_opts;
    const size_t m_max_exlvl;
    buffered::Numbers<defs::wf_t, 1> m_av_ref;
    std::vector<RefExcitsOneExlvl> m_ref_excits;
    /**
     * work space for computing connections between reference and contributing ONVs
     */
    conn::FrmOnv m_conn;

    RefExcits(const fciqmc_config::RefExcits& opts, size_t nsite) :
            Archivable("ref_excits", opts.m_archivable),
        m_opts(opts), m_max_exlvl(opts.m_max_exlvl), m_av_ref({1}), m_conn(nsite) {
        m_ref_excits.reserve(m_max_exlvl + 1);
        for (size_t i=1ul; i<=m_max_exlvl; ++i) m_ref_excits.emplace_back(i, i, 1);
    }

    void make_contribs(const fields::mbf_t& mbf, const fields::mbf_t& ref_mbf, const defs::wf_t& contrib, const size_t& ipart) {
        m_conn.connect(ref_mbf, mbf);
        auto nop = m_conn.size();
        if (!nop) m_av_ref[ipart] += contrib;
        else if (nop<=m_max_exlvl) m_ref_excits[nop - 1].make_contribs(m_conn, contrib, ipart);
    }

    bool all_stores_empty() const {
        for (const auto& it: m_ref_excits) if(it.m_hwm) return false;
        return true;
    };

private:
    /*
     *
    void save(hdf5::GroupWriter& parent) const {
    }
     */
    void load_fn(hdf5::GroupReader &parent) override {

    }

    void save_fn(hdf5::GroupWriter &parent) override {
        hdf5::GroupWriter gw("ref_excits", parent);
        defs::wf_t av_ref;
        av_ref = mpi::all_sum(m_av_ref[0]);
        gw.save("0", av_ref);
        for (auto& it: m_ref_excits) it.save(gw);
    }

};

#endif //M7_REFEXCITS_H
