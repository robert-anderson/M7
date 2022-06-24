//
// Created by Robert J. Anderson on 23/07/2021.
//

#include "FrmOnvConnection.h"
#include "M7_lib/util/Exsig.h"

FrmOnvConnection::FrmOnvConnection(const sys::frm::Size& sites):
        m_ann(sites), m_cre(sites),
        m_ndataword(integer::divceil(sites.m_nspinorb, Buffer::c_nbit_word)),
        m_dataword_phases(m_ndataword){
    if (m_ndataword) m_dataword_phases[0] = false;
}

FrmOnvConnection::FrmOnvConnection(sys::Size size) : FrmOnvConnection(size.m_frm) {
    size.require_pure_frm();
}

FrmOnvConnection::FrmOnvConnection(const FrmOnvField& mbf) : FrmOnvConnection(mbf.m_basis){}

void FrmOnvConnection::connect(const FrmOnvField& src, const FrmOnvField &dst) {
    DEBUG_ASSERT_TRUE(src.m_basis==dst.m_basis, "src and dst ONVs are incompatible");
    DEBUG_ASSERT_FALSE(src.is_zero(), "should not be computing connection from zero ONV");
    DEBUG_ASSERT_FALSE(dst.is_zero(), "should not be computing connection to zero ONV");
    clear();

    uint_t src_work, dst_work, work;
    for (uint_t idataword = 0ul; idataword < src.m_dsize; ++idataword) {
        const auto bit_offset = idataword * Buffer::c_nbit_word;
        src_work = src.get_dataword(idataword);
        dst_work = dst.get_dataword(idataword);
        work = src_work & ~dst_work;
        while (work) m_ann.add(bit::next_setbit(work) + bit_offset);
        work = dst_work & ~src_work;
        while (work) m_cre.add(bit::next_setbit(work) + bit_offset);
    }
    DEBUG_ASSERT_TRUE(m_cre.is_valid(), "creation operators are not unique and in ascending order");
    DEBUG_ASSERT_TRUE(m_ann.is_valid(), "annihilation operators are not unique and in ascending order");
}

bool FrmOnvConnection::connect(const FrmOnvField &src, const FrmOnvField &dst, FrmOps &com) {
    DEBUG_ASSERT_TRUE(src.m_basis==dst.m_basis, "src and dst ONVs are incompatible");
    connect(src, dst);
    com.clear();
    uint_t nperm = 0ul;

    auto ann_iter = m_ann.cbegin();
    auto cre_iter = m_cre.cbegin();
    const auto ann_end = m_ann.cend();
    const auto cre_end = m_cre.cend();

    uint_t src_work, dst_work, work;
    for (uint_t idataword = 0ul; idataword < src.m_dsize; ++idataword) {
        src_work = src.get_dataword(idataword);
        dst_work = dst.get_dataword(idataword);
        work = src_work & dst_work;
        while (work) {
            auto setbit = bit::next_setbit(work) + idataword * Buffer::c_nbit_word;
            while (ann_iter != ann_end && *ann_iter < setbit) {
                // an annihilation operator has been passed in the iteration over common indices
                ann_iter++;
                nperm += com.size();
            }
            while (cre_iter != cre_end && *cre_iter < setbit) {
                // a creation operator has been passed in the iteration over common indices
                cre_iter++;
                nperm += com.size();
            }
            com.add(setbit);
        }
    }
    while (ann_iter != ann_end) {
        ann_iter++;
        nperm += com.size();
    }
    while (cre_iter != cre_end) {
        cre_iter++;
        nperm += com.size();
    }
    DEBUG_ASSERT_TRUE(m_cre.is_valid(), "creation operators are not unique and in ascending order");
    DEBUG_ASSERT_TRUE(m_ann.is_valid(), "annihilation operators are not unique and in ascending order");
    return nperm & 1ul;
}


void FrmOnvConnection::apply(const FrmOnvField &src, FrmOnvField &dst) const {
    DEBUG_ASSERT_TRUE(src.m_basis==dst.m_basis, "src and dst ONVs are incompatible");
    DEBUG_ASSERT_FALSE(src.is_zero(), "should not be computing connection from zero ONV");
    DEBUG_ASSERT_TRUE(m_cre.is_valid(), "creation operators are not unique and in ascending order");
    DEBUG_ASSERT_TRUE(m_ann.is_valid(), "annihilation operators are not unique and in ascending order");
    const auto nann = m_ann.size();
    const auto ncre = m_cre.size();
    DEBUG_ASSERT_TRUE(m_ann.all_occ(src), "not all annihilation indices are occupied in src ONV");
    DEBUG_ASSERT_TRUE(m_cre.all_vac(src), "not all creation indices are vacant in src ONV");
    dst = src;
    for (uint_t i = 0ul; i < nann; ++i) dst.clr(m_ann[i]);
    for (uint_t i = 0ul; i < ncre; ++i) dst.set(m_cre[i]);
    DEBUG_ASSERT_EQ(src.nsetbit(), dst.nsetbit(),
                    "currently, all excitations are particle number conserving, so something has gone wrong here");
    DEBUG_ASSERT_TRUE(m_ann.all_vac(dst), "not all annihilation indices are vacant in dst ONV");
    DEBUG_ASSERT_TRUE(m_cre.all_occ(dst), "not all creation indices are occupied in dst ONV");
}

bool FrmOnvConnection::apply(const FrmOnvField &src, FrmOps &com) const {
    DEBUG_ASSERT_TRUE(m_cre.is_valid(), "creation operators are not unique and in ascending order");
    DEBUG_ASSERT_TRUE(m_ann.is_valid(), "annihilation operators are not unique and in ascending order");
    com.clear();
    uint_t nperm = 0ul;

    auto ann_iter = m_ann.cbegin();
    auto cre_iter = m_cre.cbegin();
    const auto ann_end = m_ann.cend();
    const auto cre_end = m_cre.cend();

    for (uint_t idataword = 0ul; idataword < src.m_dsize; ++idataword) {
        auto work = src.get_dataword(idataword);
        while (work) {
            auto setbit = bit::next_setbit(work) + idataword * Buffer::c_nbit_word;
            if (ann_iter != ann_end && setbit == *ann_iter) {
                ann_iter++;
                nperm += com.size();
                continue;
            }
            DEBUG_ASSERT_TRUE((cre_iter == cre_end) || (setbit != *cre_iter),
                              "can't create fermion in an occupied orbital");
            while (cre_iter != cre_end && *cre_iter < setbit) {
                cre_iter++;
                nperm += com.size();
            }
            com.add(setbit);
        }
    }
    while (cre_iter != cre_end) {
        cre_iter++;
        nperm += com.size();
    }
    DEBUG_ASSERT_TRUE(m_cre.is_valid(), "creation operators are not unique and in ascending order");
    DEBUG_ASSERT_TRUE(m_ann.is_valid(), "annihilation operators are not unique and in ascending order");
    return nperm & 1ul;
}

bool FrmOnvConnection::apply(const FrmOnvField &src, FrmOnvField &dst, FrmOps &com) const {
    apply(src, dst);
    return apply(src, com);
}

void FrmOnvConnection::clear() {
    m_cre.clear();
    m_ann.clear();
}

const uintv_t &FrmOnvConnection::ann() const {
    return m_ann.inds();
}

const uintv_t &FrmOnvConnection::cre() const {
    return m_cre.inds();
}

void FrmOnvConnection::update_dataword_phases(const FrmOnvField &src) const {
    for (uint_t idataword = 1ul; idataword < m_ndataword; ++idataword) {
        auto prev_dataword = src.get_dataword(idataword - 1);
        bool phase = bit::nsetbit(prev_dataword) & 1ul;
        m_dataword_phases[idataword] = (m_dataword_phases[idataword - 1] != phase);
    }
}

bool FrmOnvConnection::independent_phase(const FrmOnvField &src, const uint_t &ibit) const {
    DEBUG_ASSERT_LT(ibit, src.m_basis.m_nspinorb, "spin orbital index OOB");
    auto idataword = ibit / Buffer::c_nbit_word;
    DEBUG_ASSERT_LT(idataword, m_ndataword, "dataword index OOB");
    auto ibit_in_word = ibit - idataword * Buffer::c_nbit_word;
    return m_dataword_phases[idataword] ^
           (bit::nsetbit_before(src.get_dataword(idataword), ibit_in_word) & 1ul);
}

bool FrmOnvConnection::phase(const FrmOnvField &src) const {
    DEBUG_ASSERT_TRUE(m_ann.all_occ(src), "not all annihilation indices are occupied in src ONV");
    DEBUG_ASSERT_TRUE(m_cre.all_vac(src), "not all creation indices are vacant in src ONV");
    DEBUG_ASSERT_TRUE(m_cre.is_valid(), "creation operators are not unique and in ascending order");
    DEBUG_ASSERT_TRUE(m_ann.is_valid(), "annihilation operators are not unique and in ascending order");
    bool out = false;
    update_dataword_phases(src);
    auto ann_iter = m_ann.cbegin();
    auto cre_iter = m_cre.cbegin();
    const auto ann_end = m_ann.cend();
    const auto cre_end = m_cre.cend();

    while (ann_iter != ann_end || cre_iter != cre_end) {
        if (cre_iter!=cre_end && (ann_iter==ann_end || *cre_iter < *ann_iter)) {
            bool ann_remain_phase = std::distance(ann_iter, ann_end)&1l;
            out ^= independent_phase(src, *cre_iter++) ^ ann_remain_phase;
        }
        else {
            DEBUG_ASSERT_FALSE(ann_iter==ann_end, "should not have reached the end of annihilation indices");
            DEBUG_ASSERT_TRUE(cre_iter==cre_end || *ann_iter < *cre_iter, "invalid connection");
            out ^= independent_phase(src, *ann_iter++);
        }
    }
    out ^= (m_ann.size()/2)&1ul;
    // m*n is odd if both are odd
    out ^= (m_ann.size()&1ul) & (m_cre.size()&1ul);
    return out;
}

uint_t FrmOnvConnection::exsig() const {
    return exsig::encode(m_cre.size(), m_ann.size(), 0ul, 0ul);
}

uint_t FrmOnvConnection::exsig(const uint_t& nop_insert) const {
    return exsig::encode(m_cre.size() + nop_insert, m_ann.size() + nop_insert, 0ul, 0ul);
}