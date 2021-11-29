//
// Created by anderson on 11/29/21.
//

#include "Mbf.h"


void mbf::set_from_def(field::FrmOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef) {
    if (def.m_frm.get().empty()) return;
    REQUIRE_LT(idef, def.m_frm.get().size(), "MBF definition index OOB");
    mbf.zero();
    auto& definds = def.m_frm.get()[idef];
    if (definds.size() == mbf.m_nspinorb) {
        // assume the determinant is specified as a bit string
        for (size_t i=0ul; i<definds.size(); ++i) {
            auto flag = definds[i];
            REQUIRE_LE(flag, 1ul, "this MBF definition is being treated as a bitstring but value is not 0 or 1");
            if (flag) mbf.set(i);
        }
    }
    else {
        // assume the determinant is specified as a vector of set spin orbital indices
        for (auto& ind: definds){
            REQUIRE_LT(ind, mbf.m_nspinorb, "spin orbital index OOB");
            mbf.set(ind);
        }
    }
}

void mbf::set_from_def(field::FrmBosOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef) {
    set_from_def(mbf.m_frm, def, idef);
    set_from_def(mbf.m_bos, def, idef);
}

void mbf::set_from_def(field::BosOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef) {
    if (def.m_bos.get().empty()) return;
    REQUIRE_LT(idef, def.m_bos.get().size(), "MBF definition index OOB");
    mbf.zero();
    auto& definds = def.m_bos.get()[idef];
    REQUIRE_EQ(definds.size(), mbf.m_nelement,
               "length of input vector does not match number of boson modes");
    size_t i = 0ul;
    for (auto &occ: definds) mbf[i++] = occ;
}
