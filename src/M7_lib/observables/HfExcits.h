//
// Created by Robert J. Anderson on 03/03/2021.
//

#ifndef M7_HFEXCITS_H
#define M7_HFEXCITS_H

#include <M7_lib/connection/Connections.h>
#include <M7_lib/conf/Conf.h>
#include <M7_lib/io/Archivable.h>
#include <M7_lib/basis/Suites.h>
#include <M7_lib/field/Fields.h>
#include <M7_lib/table/BufferedTable.h>
#include <M7_lib/table/BufferedFields.h>
#include <M7_lib/mae/MaeTable.h>

/**
 * Averages weights on the MBFs local to the HF MBF in the Hilbert space.
 * This MBF need not be a true HF determinant with a Brillouin theorem, so it is sometimes referred to in the code as
 * HF-like
 */
struct HfExcitsOneExsig : buffered::MappedTable<MaeRow> {
protected:
    /**
     * work space for converting between indexing vector type and the stored key type of the MEV tables
     */
    mutable buffered::MaeInds m_working_inds;

public:
    HfExcitsOneExsig(OpSig exsig, uint_t nroot);

    MaeRow& lookup(const conn::FrmOnv& key);

    MaeRow& insert(const conn::FrmOnv& key);

    strv_t h5_field_names() const;

    using Table<MaeRow>::save;
    void save(const hdf5::NodeWriter& gw) const;

    void make_contribs(const conn::FrmOnv& conn, const wf_t& contrib, uint_t iroot);
};


struct HfExcits : Archivable {
    /**
     * relevant section of configuration document
     */
    const conf::HfExcits& m_opts;
    /**
     * averaged weights on HF-like MBF
     */
    buffered::Numbers<wf_t, 1> m_av_hf;
    /**
     * all one-exsig objects
     */
    std::array<std::unique_ptr<HfExcitsOneExsig>, opsig::c_ndistinct> m_hf_excits;
    /**
     * excitation signatures for which HF excits are currently accumulated
     */
    v_t<OpSig> m_active_exsigs;
    /**
     * work space for computing connections between reference and contributing ONVs
     */
    conn::Mbf m_conn;

    HfExcits(const conf::HfExcits& opts, sys::Size extents, uint_t nroot);

    void make_contribs(const conn::FrmOnv& conn, const wf_t& contrib, uint_t iroot);

    void make_contribs(const conn::FrmBosOnv& conn, const wf_t& contrib, uint_t iroot);

    void make_contribs(const conn::BosOnv& /*conn*/, const wf_t& /*contrib*/, uint_t /*iroot*/){
        ABORT("not yet implemented");
    }

    void make_contribs(const field::Mbf& mbf, const field::Mbf& hf_mbf, const wf_t& contrib, uint_t iroot);

    bool all_stores_empty() const;

    operator bool() const;

private:

    void load_fn(const hdf5::NodeReader& /*parent*/) override {}

    void save_fn(const hdf5::NodeWriter& parent) override;
};

#endif //M7_HFEXCITS_H
