//
// Created by rja on 08/08/2021.
//

#include "FrmOnvField.h"


FrmOnvField::FrmOnvField(Row *row, BasisData bd, std::string name) :
        base_t(row, {{2, bd.m_nsite},
                     {"spin channel", "site"}}, name),
        m_bd(bd), m_nsite(bd.m_nsite), m_nspinorb(m_format.m_nelement),
        m_decoded(*this, bd.m_frm_abgrp_map),
        m_dsize_spin_channel(integer_utils::divceil(m_nsite, defs::nbit_word)),
        m_nbit_in_last_alpha_dataword(m_nsite%defs::nbit_word){
    bd.require_pure_frm();
}

FrmOnvField::FrmOnvField(Row *row, size_t nsite, std::string name) :
    FrmOnvField(row, {nsite, 0ul}, name){}

FrmOnvField::FrmOnvField(const FrmOnvField &other) :
        base_t(other), m_bd(other.m_bd), m_nsite(other.m_nsite), m_nspinorb(other.m_nspinorb),
        m_decoded(*this, other.m_bd.m_frm_abgrp_map),
        m_dsize_spin_channel(other.m_dsize_spin_channel),
        m_nbit_in_last_alpha_dataword(other.m_nbit_in_last_alpha_dataword){}

FrmOnvField &FrmOnvField::operator=(std::pair<const defs::inds &, const defs::inds &> setbits) {
    // prezero the element
    zero();
    for (const auto &ind: setbits.first) set(ind);
    for (const auto &ind: setbits.second) set(ind+m_nsite);
    return *this;
}

bool FrmOnvField::operator==(const defs::inds &inds) const {
    if (nsetbit()!=inds.size()) return false;
    for (const auto& ind: inds) if (!get(ind)) return false;
    return true;
}

int FrmOnvField::ms2() const {
    int spin = 0;
    size_t work;
    for (size_t idataword = 0; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            if (ibit < m_nsite) ++spin;
            else if (ibit >= nbit()) return spin;
            else --spin;
        }
    }
    return spin;
}

int FrmOnvField::nalpha() const {
    // number of electrons occupying spinors in the alpha spin channel
    int nalpha = 0;
    size_t work;
    for (size_t idataword = 0ul; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            if (ibit >= m_nsite) return nalpha;
            nalpha++;
        }
    }
    return nalpha;
}

size_t FrmOnvField::site_nocc(const size_t &isite) const {
    return get({0, isite}) + get({1, isite});
}

bool FrmOnvField::all_sites_single_occ() const {
    for (size_t isite = 0ul; isite<m_nsite; ++isite) if (site_nocc(isite)!=1ul) return false;
    return true;
}

std::string FrmOnvField::to_string() const {
    std::string res;
    res += "(";
    res.reserve(nbit() + 3);
    size_t i = 0ul;
    for (; i < m_nsite; ++i)
        res += get(i) ? "1" : "0";
    res += ","; // spin channel delimiter
    for (; i < nbit(); ++i)
        res += get(i) ? "1" : "0";
    res += ")";
    return res;
}

const size_t &FrmOnvField::nsite() const {
    return m_nsite;
}

size_t FrmOnvField::isite(const size_t &ibit, const size_t &nsite) {
    return ibit < nsite ? ibit : ibit - nsite;
}

size_t FrmOnvField::ispin(const size_t &ibit, const size_t &nsite) {
    return ibit >= nsite;
}

size_t FrmOnvField::ibit(const size_t &ispin, const size_t &isite, const size_t &nsite) {
    return ispin ? nsite+isite : isite;
}

size_t FrmOnvField::ibit(std::pair<size_t, size_t> pair, const size_t &nsite) {
    return ibit(pair.first, pair.second, nsite);
}


size_t FrmOnvField::isite(const size_t &ibit) const {
    return isite(ibit, m_nsite);
}

size_t FrmOnvField::ispin(const size_t &ibit) const {
    return ispin(ibit, m_nsite);
}

size_t FrmOnvField::ibit(const size_t &ispin, const size_t &isite) const {
    return ibit(ispin, isite, m_nsite);
}

size_t FrmOnvField::ibit(std::pair<size_t, size_t> &pair) const {
    return ibit(pair.first, pair.second);
}

void FrmOnvField::set(const size_t &bit_offset, const defs::inds &setbits) {
    for(auto& i: setbits) set(bit_offset+i);
}

void FrmOnvField::set(const size_t &site_offset, const defs::inds &setbits_alpha, const defs::inds &setbits_beta) {
    for(auto& i: setbits_alpha) set({0, site_offset+i});
    for(auto& i: setbits_beta) set({1, site_offset+i});
}

void FrmOnvField::set(const defs::inds &setbits_alpha, const defs::inds &setbits_beta) {
    zero();
    for(auto& i: setbits_alpha) set({0, i});
    for(auto& i: setbits_beta) set({1, i});
}

void FrmOnvField::set_spins(const defs::inds &alpha_sites) {
    DEBUG_ASSERT_LE(alpha_sites.size(), m_nsite, "can't have more spins than sites");
    zero();
    auto it = alpha_sites.cbegin();
    for (size_t isite=0ul; isite<m_nsite; ++isite){
        if (it==alpha_sites.cend() || *it>isite) set({1, isite});
        else {
            set({0, isite});
            ++it;
        };
    }
}

void FrmOnvField::put_spin_channel(const size_t &ispin, bool set) {
    auto ibegin = ibit(ispin, 0);
    put_range(ibegin, ibegin+m_nsite, set);
}

void FrmOnvField::clr(const size_t &bit_offset, const defs::inds &clrbits) {
    for(auto& i: clrbits) clr(bit_offset+i);
}

void FrmOnvField::clr(const size_t &site_offset, const defs::inds &clrbits_alpha, const defs::inds &clrbits_beta) {
    for(auto& i: clrbits_alpha) clr({0, site_offset+i});
    for(auto& i: clrbits_beta) clr({1, site_offset+i});
}

size_t FrmOnvField::get_alpha_dataword(size_t idataword) const {
    DEBUG_ASSERT_LT(idataword, m_dsize_spin_channel, "dataword index OOB");
    size_t *dptr = reinterpret_cast<size_t *>(begin());
    auto tmp = dptr[idataword];
    if (idataword + 1 == m_dsize_spin_channel) {
        tmp = bit_utils::truncate(tmp, m_nbit_in_last_alpha_dataword);
    }
    return tmp;
}

size_t FrmOnvField::get_beta_dataword(size_t idataword) const {
    /*
     * suppose we are storing an nsite=19 FrmOnv with 8-bit dataword:
     *
     * |aaaaaaaa|aaaaaaaa|aaabbbbb|bbbbbbbb|bbbbbb--|
     *
     * shift to the left
     *                   |bbbbb---|
     * shift to the right
     *                            |-----bbb|
     * then combine
     */
    DEBUG_ASSERT_LT(idataword, m_dsize_spin_channel, "dataword index OOB");
    size_t *dptr = reinterpret_cast<size_t *>(begin());
    auto left = dptr[m_dsize_spin_channel-1+idataword];
    left >>= m_nbit_in_last_alpha_dataword;
    auto nbit_in_left = defs::nbit_word-m_nbit_in_last_alpha_dataword;
    bit_utils::truncate(left, m_nbit_in_last_alpha_dataword);
    auto right = dptr[m_dsize_spin_channel+idataword];
    right <<= nbit_in_left;
    return left|right;
}

size_t FrmOnvField::nopenshell() const {
    size_t count = 0ul;
    auto fn = [&count](size_t isite){++count;};
    foreach_openshell(fn);
    return count;
}

size_t FrmOnvField::nopenshell_alpha() const {
    size_t count = 0ul;
    auto fn = [&count](size_t isite){++count;};
    foreach_openshell_alpha(fn);
    return count;
}
