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

struct RefExcitsOneExsig : BufferedTable<MaeRow, true> {
    /**
     * work space for converting between indexing vector type and the stored key type of the MEV tables
     */
    buffered::MaeInds m_working_inds;

    RefExcitsOneExsig(size_t exsig, size_t nvalue, size_t nbucket = 100);

    LookupResult operator[](const conn::FrmOnv &key);

    size_t insert(const conn::FrmOnv &key);

    std::vector<std::string> h5_field_names() const;

    using Table<MaeRow>::save;
    void save(hdf5::GroupWriter& gw) const;

    void make_contribs(const conn::FrmOnv& conn, const defs::wf_t& contrib, const size_t& ipart);
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

    RefExcits(const fciqmc_config::RefExcits& opts, size_t nsite);

    void make_contribs(const conn::FrmOnv& conn, const defs::wf_t& contrib, const size_t& ipart);

    void make_contribs(const conn::FrmBosOnv& conn, const defs::wf_t& contrib, const size_t& ipart);

    void make_contribs(const field::Mbf& mbf, const field::Mbf& ref_mbf, const defs::wf_t& contrib, const size_t& ipart);

    bool all_stores_empty() const;

    operator bool() const;

private:

    void load_fn(hdf5::GroupReader &parent) override;

    void save_fn(hdf5::GroupWriter &parent) override;
};

#endif //M7_REFEXCITS_H
