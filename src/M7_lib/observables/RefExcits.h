//
// Created by rja on 03/03/2021.
//

#ifndef M7_REFEXCITS_H
#define M7_REFEXCITS_H

#include <M7_lib/connection/Connections.h>
#include <M7_lib/config/FciqmcConfig.h>
#include <M7_lib/io/Archivable.h>
#include <M7_lib/basis/Suites.h>
#include <M7_lib/field/Fields.h>
#include <M7_lib/table/BufferedTable.h>
#include <M7_lib/table/BufferedFields.h>
#include <M7_lib/mae/MaeTable.h>

struct RefExcitsOneExsig : BufferedTable<MaeRow, true> {
    /**
     * work space for converting between indexing vector type and the stored key type of the MEV tables
     */
    buffered::MaeInds m_working_inds;

    RefExcitsOneExsig(size_t exsig, size_t nroot, size_t nbucket = 100);

    LookupResult operator[](const conn::FrmOnv &key);

    size_t insert(const conn::FrmOnv &key);

    std::vector<std::string> h5_field_names() const;

    using Table<MaeRow>::save;
    void save(hdf5::GroupWriter& gw) const;

    void make_contribs(const conn::FrmOnv& conn, const defs::wf_t& contrib, const size_t& iroot);
};


struct RefExcits : Archivable {
    const conf::RefExcits& m_opts;
    buffered::Numbers<defs::wf_t, 1> m_av_ref;
    std::array<std::unique_ptr<RefExcitsOneExsig>, defs::nexsig> m_ref_excits;
    defs::inds m_active_exsigs;
    /**
     * work space for computing connections between reference and contributing ONVs
     */
    conn::Mbf m_conn;

    RefExcits(const conf::RefExcits& opts, sys::Size extents, size_t nroot);

    void make_contribs(const conn::FrmOnv& conn, const defs::wf_t& contrib, const size_t& iroot);

    void make_contribs(const conn::FrmBosOnv& conn, const defs::wf_t& contrib, const size_t& iroot);

    void make_contribs(const conn::BosOnv& conn, const defs::wf_t& contrib, const size_t& iroot);

    void make_contribs(const field::Mbf& mbf, const field::Mbf& ref_mbf, const defs::wf_t& contrib, const size_t& iroot);

    bool all_stores_empty() const;

    operator bool() const;

private:

    void load_fn(hdf5::GroupReader &parent) override;

    void save_fn(hdf5::GroupWriter &parent) override;
};

#endif //M7_REFEXCITS_H
