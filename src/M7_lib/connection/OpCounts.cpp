//
// Created by rja on 26/03/23.
//

#include "OpCounts.h"

void OpCounts::set(const field::FrmOnv& src, const field::FrmOnv& dst) {
    uint_t src_work, dst_work, work;
    m_nfrm_cre = 0ul;
    m_nfrm_ann = 0ul;
    for (uint_t idataword = 0ul; idataword < src.m_dsize; ++idataword) {
        src_work = src.get_dataword(idataword);
        dst_work = dst.get_dataword(idataword);
        work = src_work & ~dst_work;
        while (work) m_nfrm_ann += (bit::next_setbit(work), 1ul);
        work = dst_work & ~src_work;
        while (work) m_nfrm_cre += (bit::next_setbit(work), 1ul);
    }
}

void OpCounts::set(const field::BosOnv& src, const field::BosOnv& dst) {
    m_nbos_cre = 0ul;
    m_nbos_ann = 0ul;
    for (uint_t imode = 0ul; imode < src.nelement(); ++imode) {
        if (src[imode] > dst[imode]) ++m_nbos_ann;
        else if (src[imode] < dst[imode]) ++m_nbos_cre;
    }
}

void OpCounts::set(const field::FrmBosOnv& src, const field::FrmBosOnv& dst) {
    set(src.m_frm, dst.m_frm);
    set(src.m_bos, dst.m_bos);
}
