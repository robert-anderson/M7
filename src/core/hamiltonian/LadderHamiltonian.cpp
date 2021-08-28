//
// Created by rja on 05/11/2020.
//

#include "LadderHamiltonian.h"

LadderHamiltonian::LadderHamiltonian(size_t nmode, size_t nboson_max, std::string fname) :
        m_nboson_max(nboson_max), m_nmode(nmode), m_v(nmode), m_v_unc(nmode, 0.0),
        m_contribs_0010(exsig_utils::ex_0010), m_contribs_0001(exsig_utils::ex_0001),
        m_contribs_1110(exsig_utils::ex_1110), m_contribs_1101(exsig_utils::ex_1101) {
    if (!m_nboson_max) return;

    defs::inds inds(3);
    defs::ham_t value;
    EbdumpFileReader file_reader(fname);
    REQUIRE_EQ_ALL(file_reader.m_nspatorb, m_nmode, "expected number of boson modes not found in file");

    log::info("Reading boson ladder coupled and uncoupled coefficients from file \"" + file_reader.m_fname + "\"...");
    while (file_reader.next(inds, value)) {
        if (consts::float_is_zero(value)) continue;
        auto ranksig = file_reader.ranksig(inds);
        auto exsig = file_reader.exsig(inds, ranksig);
        if (ranksig==exsig_utils::ex_0010){
            m_contribs_0010.set_nonzero(exsig);
            m_contribs_0001.set_nonzero(exsig_utils::hermconj(exsig));
            m_v_unc[inds[0]] = value;
        }
        else {
            DEBUG_ASSERT_EQ(ranksig, exsig_utils::ex_1110, "ranksig should be either 0010 or 1110");
            m_contribs_1110.set_nonzero(exsig);
            m_contribs_1101.set_nonzero(exsig_utils::hermconj(exsig));
            m_v.set(inds[0], inds[1], inds[2], value);
        }
    }
    log_data();
}

defs::ham_t LadderHamiltonian::get_element(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    if (conn.m_bos.size() != 1ul) return 0.0;
    if (!conn.m_frm.kramers_conserve()) return 0.0;
    bool cre = conn.m_bos.m_cre.size();

    const auto imode = cre ? conn.m_bos.m_cre[0].m_imode : conn.m_bos.m_ann[0].m_imode;
    const auto com = size_t(onv.m_bos[imode]) - (cre ? 0 : 1);

    const auto occ_fac = std::sqrt(com + 1);
    switch (conn.m_frm.size()) {
        case 0: {
            defs::ham_t res = m_v_unc[imode];
            // fermion ONVs do not differ, so sum over occupied spin orbitals
            auto fn = [&](const size_t &ibit) {
                auto isite = ibit % m_nmode;
                res += m_v.get(imode, isite, isite);
            };
            onv.m_frm.foreach(fn);
            return res * occ_fac;
        }
        case 2: {
            DEBUG_ASSERT_TRUE(onv.m_frm.get(conn.m_frm.m_ann[0]), "annihilated op not occupied in ONV")
            DEBUG_ASSERT_FALSE(onv.m_frm.get(conn.m_frm.m_cre[0]), "created op occupied in ONV")
            auto isite = conn.m_frm.m_cre[0] % m_nmode;
            auto jsite = conn.m_frm.m_ann[0] % m_nmode;
            /*
             * respect hermitian conjugation of the fermion-boson operator product: 1110 (boson creation) is the
             * conventionally non-conjugated term
             */
            auto element = cre ? m_v.get(imode, isite, jsite) : m_v.get(imode, jsite, isite);
            element*=occ_fac;
            return conn.m_frm.phase(onv.m_frm) ? -element : element;
        }
        default:
            return 0;
    }
}

void LadderHamiltonian::log_data() const {
    if (!m_contribs_0010.is_nonzero(exsig_utils::ex_0010))
        log::info("0010 uncoupled boson ladder hamiltonian term has no contributions");
    if (!m_contribs_0001.is_nonzero(exsig_utils::ex_0001))
        log::info("0001 uncoupled boson ladder hamiltonian term has no contributions");

    if (!m_contribs_1110.is_nonzero(exsig_utils::ex_0010))
        log::info("1110 fermion-coupled boson ladder term has no 0010 contributions");
    if (!m_contribs_1110.is_nonzero(exsig_utils::ex_1110))
        log::info("1110 fermion-coupled boson ladder term has no 1110 contributions");
    if (!m_contribs_1101.is_nonzero(exsig_utils::ex_0001))
        log::info("1101 fermion-coupled boson ladder term has no 0001 contributions");
    if (!m_contribs_1101.is_nonzero(exsig_utils::ex_1101))
        log::info("1101 fermion-coupled boson ladder term has no 1101 contributions");
}
