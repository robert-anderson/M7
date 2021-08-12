//
// Created by rja on 03/03/2021.
//

#ifndef M7_REFEXCITS_H
#define M7_REFEXCITS_H

#include <src/core/connection/Connections.h>
#include <src/core/config/FciqmcConfig.h>
#include <src/core/io/Archivable.h>
#include <src/core/basis/Suites.h>
#include "src/core/field/Fields.h"
#include "src/core/table/BufferedTable.h"
#include "src/core/table/BufferedFields.h"
#include "src/core/mae/MaeTable.h"

struct RefExcitsOneExsig : BufferedTable<MaeRow<defs::wf_t>, true> {
    /**
     * work space for converting between indexing vector type and the stored key type of the MEV tables
     */
    buffered::MaeInds m_working_inds;

    RefExcitsOneExsig(size_t exsig, size_t nvalue, size_t nbucket = 100) :
            BufferedTable<MaeRow<defs::wf_t>, true>("average coefficients", {{exsig, nvalue}, nbucket}),
            m_working_inds(exsig) {
    }

    LookupResult operator[](const conn::FrmOnv &key) {
        m_working_inds = key;
        return MappedTable<MaeRow<defs::wf_t>>::operator[](m_working_inds);
    }

    size_t insert(const conn::FrmOnv &key) {
        m_working_inds = key;
        return MappedTable<MaeRow<defs::wf_t>>::insert(m_working_inds);
    }

    std::vector<std::string> h5_field_names() const {
        return {m_row.m_inds.m_name, m_row.m_values.m_name};
    }

    using Table<MaeRow<defs::wf_t>>::save;
    void save(hdf5::GroupWriter& gw) const {
        Table<MaeRow<defs::wf_t>>::save(gw, m_working_inds.get_exsig_string(), h5_field_names());
    }

    void make_contribs(const conn::FrmOnv& conn, const defs::wf_t& contrib, const size_t& ipart) {
        DEBUG_ASSERT_EQ(conn.exsig(), m_working_inds.m_exsig, "incompatible connection");
        auto irow = *(*this)[conn];
        if (irow==~0ul) irow = insert(conn);
        m_row.jump(irow);
        m_row.m_values[ipart]+=contrib;
    }
};


struct RefExcits : Archivable {
    const fciqmc_config::RefExcits& m_opts;
    buffered::Numbers<defs::wf_t, 1> m_av_ref;
    std::array<std::unique_ptr<RefExcitsOneExsig>, defs::nexsig> m_ref_excits;
    defs::inds m_active_exsigs;
    /**
     * work space for computing connections between reference and contributing ONVs
     */
    conn::Mbf m_conn;

    RefExcits(const fciqmc_config::RefExcits& opts, size_t nsite) :
            Archivable("ref_excits", opts.m_archivable),
        m_opts(opts), m_av_ref({1}), m_conn(nsite) {
        for (size_t iexlvl=1ul; iexlvl<opts.m_max_exlvl; ++iexlvl){
            auto exsig = conn_utils::exsig(iexlvl, iexlvl, 0, 0);
            m_active_exsigs.push_back(exsig);
            m_ref_excits[exsig] = mem_utils::make_unique<RefExcitsOneExsig>(exsig, 1);
        }
    }

private:
    void make_contribs(const conn::FrmOnv& conn, const defs::wf_t& contrib, const size_t& ipart) {
        auto exsig = conn.exsig();
        if (!exsig) m_av_ref[ipart] += contrib;
        if (exsig==~0ul || !m_ref_excits[exsig]) return;
        m_ref_excits[exsig]->make_contribs(conn, contrib, ipart);
    }

    void make_contribs(const conn::FrmBosOnv& conn, const defs::wf_t& contrib, const size_t& ipart) {
        make_contribs(conn.m_frm, contrib, ipart);
    }

public:

    void make_contribs(const field::Mbf& mbf, const field::Mbf& ref_mbf, const defs::wf_t& contrib, const size_t& ipart) {
        m_conn.connect(ref_mbf, mbf);
        make_contribs(m_conn, contrib, ipart);
    }

    bool all_stores_empty() const {
        for (const auto& i: m_active_exsigs) {
            DEBUG_ASSERT_TRUE(m_ref_excits[i].get(), "active exsig was not allocated!");
            if (m_ref_excits[i]->m_hwm) return false;
        }
        return true;
    };

private:

    void load_fn(hdf5::GroupReader &parent) override {

    }

    void save_fn(hdf5::GroupWriter &parent) override {
        hdf5::GroupWriter gw("ref_excits", parent);
        defs::wf_t av_ref;
        av_ref = mpi::all_sum(m_av_ref[0]);
        gw.save("0000", av_ref);
        for (const auto& i: m_active_exsigs) {
            DEBUG_ASSERT_TRUE(m_ref_excits[i].get(), "active exsig was not allocated!");
            m_ref_excits[i]->save(gw);
        }
    }
};

#endif //M7_REFEXCITS_H
