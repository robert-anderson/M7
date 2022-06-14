//
// Created by Robert J. Anderson on 08/08/2021.
//

#include "FrmOnvField.h"

FrmOnvField::FrmOnvField(Row *row, const sys::frm::Basis& basis, std::string name) :
        base_t(row, {{2, basis.m_nsite},{"spin channel", "site"}}, name),
        m_basis(basis), m_decoded(*this),
        m_dsize_spin_channel(integer_utils::divceil(size_t(m_basis.m_nsite), defs::nbit_word)){}

FrmOnvField::FrmOnvField(Row *row, const sys::Basis &basis, std::string name) :
    FrmOnvField(row, basis.m_frm, std::move(name)){
    basis.require_pure_frm();
}

FrmOnvField::FrmOnvField(Row *row, const sys::frm::Sector& sector, std::string name) :
    FrmOnvField(row, sector.m_basis, name) {}

FrmOnvField::FrmOnvField(Row *row, const sys::Sector& sector, std::string name) :
    FrmOnvField(row, sector.m_frm, name){}

FrmOnvField::FrmOnvField(const FrmOnvField &other) :
        base_t(other), m_basis(other.m_basis), m_decoded(*this),
        m_dsize_spin_channel(other.m_dsize_spin_channel){}

FrmOnvField &FrmOnvField::operator=(std::pair<const defs::inds &, const defs::inds &> setbits) {
    // prezero the element
    zero();
    for (const auto &ind: setbits.first) set(ind);
    for (const auto &ind: setbits.second) set(ind+m_basis.m_nsite);
    return *this;
}

bool FrmOnvField::operator==(const defs::inds &inds) const {
    if (nsetbit()!=inds.size()) return false;
    for (const auto& ind: inds) if (!get(ind)) return false;
    return true;
}

int FrmOnvField::ms2() const {
    int ms2 = 0;
    size_t work;
    for (size_t idataword = 0; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            if (ibit >= nbit()) return ms2;
            ms2 += m_basis.ms2(ibit);
        }
    }
    return ms2;
}

size_t FrmOnvField::site_nocc(const size_t &isite) const {
    return get({0, isite}) + get({1, isite});
}

bool FrmOnvField::all_sites_single_occ() const {
    return nopen_shell()==m_basis.m_nsite;
}

std::string FrmOnvField::to_string() const {
    std::string res;
    res += "(";
    res.reserve(nbit() + 3);
    size_t i = 0ul;
    for (; i < m_basis.m_nsite; ++i)
        res += get(i) ? "1" : "0";
    res += ","; // spin channel delimiter
    for (; i < nbit(); ++i)
        res += get(i) ? "1" : "0";
    res += ")";
    return res;
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
    DEBUG_ASSERT_LE(alpha_sites.size(), m_basis.m_nsite, "can't have more spins than sites");
    zero();
    auto it = alpha_sites.cbegin();
    for (size_t isite=0ul; isite<m_basis.m_nsite; ++isite){
        if (it==alpha_sites.cend() || *it>isite) set({1, isite});
        else {
            set({0, isite});
            ++it;
        };
    }
}

void FrmOnvField::put_spin_channel(const size_t &ispin, bool set) {
    auto ibegin = m_basis.ispinorb(ispin, 0);
    put_range(ibegin, ibegin+m_basis.m_nsite, set);
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
    const auto nsite = m_format.m_shape[1];
    auto dptr = reinterpret_cast<size_t *>(begin());
    auto tmp = dptr[idataword];
    if (idataword + 1 == m_dsize_spin_channel) {
        auto n = nsite-(idataword)*defs::nbit_word;
        tmp = bit_utils::truncate(tmp, n);
    }
    return tmp;
}

size_t FrmOnvField::get_beta_dataword(size_t idataword) const {
    DEBUG_ASSERT_LT(idataword, m_dsize_spin_channel, "dataword index OOB");
    const auto nsite = m_format.m_shape[1];
    const auto ibit_begin = nsite+idataword*defs::nbit_word;
    const auto ibit_end = std::min(2*nsite, ibit_begin+defs::nbit_word);
    const auto iword_begin = ibit_begin/defs::nbit_word;
    const auto ibit_in_word_begin = ibit_begin-iword_begin*defs::nbit_word;
    const auto iword_end = ibit_end/defs::nbit_word;
    const auto ibit_in_word_end = ibit_end-iword_end*defs::nbit_word;
    auto dptr = reinterpret_cast<size_t *>(begin());
    if (iword_begin==iword_end){
        auto mask = bit_utils::make_range_mask<size_t>(ibit_in_word_begin, ibit_in_word_end);
        return (dptr[iword_begin] & mask) >> ibit_in_word_begin;
    }
    else {
        auto mask1 = bit_utils::make_range_mask<size_t>(ibit_in_word_begin, defs::nbit_word);
        auto mask2 = bit_utils::make_range_mask<size_t>(0, ibit_in_word_end);
        auto shift1 = ibit_in_word_begin;
        auto shift2 = defs::nbit_word-ibit_in_word_begin;
        return ((dptr[iword_begin] & mask1) >> shift1) | ((dptr[iword_end] & mask2) << shift2);
    }
}

size_t FrmOnvField::nalpha() const {
    size_t count = 0ul;
    auto fn = [&count](size_t isite){++count;};
    foreach_alpha(fn);
    return count;
}

size_t FrmOnvField::nbeta() const {
    size_t count = 0ul;
    auto fn = [&count](size_t isite){++count;};
    foreach_beta(fn);
    return count;
}

size_t FrmOnvField::nopen_shell() const {
    size_t count = 0ul;
    auto fn = [&count](size_t isite){++count;};
    foreach_open_shell(fn);
    return count;
}

size_t FrmOnvField::nopen_shell_alpha() const {
    size_t count = 0ul;
    auto fn = [&count](size_t isite){++count;};
    foreach_open_shell_alpha(fn);
    return count;
}

size_t FrmOnvField::nopen_shell_beta() const {
    size_t count = 0ul;
    auto fn = [&count](size_t isite){++count;};
    foreach_open_shell_beta(fn);
    return count;
}

size_t FrmOnvField::nclosed_shell() const {
    size_t count = 0ul;
    auto fn = [&count](size_t isite){++count;};
    foreach_closed_shell(fn);
    return count;
}

size_t FrmOnvField::noccupied_site() const {
    size_t count = 0ul;
    auto fn = [&count](size_t isite){++count;};
    foreach_occupied_site(fn);
    return count;
}

size_t FrmOnvField::nunoccupied_site() const {
    return m_basis.m_nsite-noccupied_site();
}

