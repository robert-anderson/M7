//
// Created by Robert J. Anderson on 27/08/2021.
//

#include "HfExcits.h"
#include "M7_lib/util/Pointer.h"

HfExcitsOneExsig::HfExcitsOneExsig(uint_t exsig, uint_t nroot) :
        buffered::MappedTable<MaeRow>(
                logging::format("average {} HF-like excitation coefficients", exsig::to_string(exsig)),
                MaeRow{exsig, nroot}),
        m_working_inds(exsig) {}

MaeRow& HfExcitsOneExsig::lookup(const conn::FrmOnv& key) {
    m_working_inds = key;
    return MappedTable<MaeRow>::lookup(m_working_inds);
}

MaeRow& HfExcitsOneExsig::insert(const conn::FrmOnv& key) {
    m_working_inds = key;
    return MappedTable<MaeRow>::insert(m_working_inds);
}

strv_t HfExcitsOneExsig::h5_field_names() const {
    return {m_row.m_inds.m_name, m_row.m_values.m_name};
}

void HfExcitsOneExsig::save(const hdf5::NodeWriter& gw) const {
    Table<MaeRow>::save(gw, exsig::to_string(m_working_inds.m_exsig), h5_field_names());
}

void HfExcitsOneExsig::make_contribs(const conn::FrmOnv& conn, const wf_t& contrib, uint_t iroot) {
    DEBUG_ASSERT_EQ(conn.exsig(), m_working_inds.m_exsig, "incompatible connection");
    if (!this->lookup(conn)) this->insert(conn);
    m_row.m_values[iroot]+=contrib;
}

HfExcits::HfExcits(const conf::HfExcits& opts, sys::Size extents, uint_t nroot) :
        Archivable("hf_excits", opts.m_archivable),
        m_opts(opts), m_av_hf({nroot}), m_conn(extents) {
    REQUIRE_EQ_ALL(nroot, 1ul, "HF excitation averaging currently only implemented for a single root");
    for (uint_t iexlvl=1ul; iexlvl<=opts.m_max_exlvl; ++iexlvl){
        auto exsig = exsig::encode(iexlvl, iexlvl, 0, 0);
        m_active_exsigs.push_back(exsig);
        m_hf_excits[exsig] = ptr::smart::make_unique<HfExcitsOneExsig>(exsig, nroot);
    }
}

void HfExcits::make_contribs(const conn::FrmOnv& conn, const wf_t& contrib, uint_t iroot) {
    auto exsig = conn.exsig();
    if (!exsig) m_av_hf[iroot] += contrib;
    if (exsig==~0ul || !m_hf_excits[exsig]) return;
    m_hf_excits[exsig]->make_contribs(conn, contrib, iroot);
}

void HfExcits::make_contribs(const conn::FrmBosOnv& conn, const wf_t& contrib, uint_t iroot) {
    make_contribs(conn.m_frm, contrib, iroot);
}

void HfExcits::make_contribs(const field::Mbf& mbf, const field::Mbf& hf_mbf, const wf_t& contrib,
                             uint_t iroot) {
    m_conn.connect(hf_mbf, mbf);
    make_contribs(m_conn, contrib, iroot);
}

bool HfExcits::all_stores_empty() const {
    for (const auto& i: m_active_exsigs) {
        DEBUG_ASSERT_TRUE(m_hf_excits[i].get(), "active encode_exsig was not allocated!");
        if (!m_hf_excits[i]->empty()) return false;
    }
    return true;
}

HfExcits::operator bool() const {
    return !m_active_exsigs.empty();
}

void HfExcits::save_fn(const hdf5::NodeWriter& parent) {
    wf_t av_hf;
    av_hf = mpi::all_sum(m_av_hf[0]);
    hdf5::GroupWriter gw(parent, "hf_excits");
    gw.save("0000", av_hf);
    for (const auto& i: m_active_exsigs) {
        DEBUG_ASSERT_TRUE(m_hf_excits[i].get(), "active exsig was not allocated!");
        m_hf_excits[i]->save(gw);
    }
}
