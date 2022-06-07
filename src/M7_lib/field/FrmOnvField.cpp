//
// Created by Robert J. Anderson on 08/08/2021.
//

#include "FrmOnvField.h"

FrmOnvField::FrmOnvField(Row *row, const sys::frm::Basis& basis, std::string name) :
        base_t(row, {{2, basis.m_nsite},{"spin channel", "site"}}, name),
        m_basis(basis), m_decoded(*this),
        m_dsize_spin_channel(integer_utils::divceil(size_t(m_basis.m_nsite), defs::nbit_word)),
        m_nbit_in_last_alpha_dataword(m_basis.m_nsite % defs::nbit_word){}

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
        m_dsize_spin_channel(other.m_dsize_spin_channel),
        m_nbit_in_last_alpha_dataword(other.m_nbit_in_last_alpha_dataword){}

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

int FrmOnvField::nalpha() const {
    // number of electrons occupying spinors in the alpha spin channel
    int nalpha = 0;
    size_t work;
    for (size_t idataword = 0ul; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            if (m_basis.ispin(ibit)) return nalpha;
            nalpha++;
        }
    }
    return nalpha;
}

size_t FrmOnvField::site_nocc(const size_t &isite) const {
    return get({0, isite}) + get({1, isite});
}

bool FrmOnvField::all_sites_single_occ() const {
    for (size_t isite = 0ul; isite<m_basis.m_nsite; ++isite) if (site_nocc(isite)!=1ul) return false;
    return true;
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
    auto dptr = reinterpret_cast<size_t *>(begin());
    auto tmp = dptr[idataword];
    if (idataword + 1 == m_dsize_spin_channel) {
        tmp = bit_utils::truncate(tmp, m_nbit_in_last_alpha_dataword);
    }
    return tmp;
}

size_t FrmOnvField::get_beta_dataword(size_t idataword) const {
    /*
     * suppose we are storing an nsite=19 FrmOnv with 8-bit dataword:
     *                       00000 00011111 111222
     * |aaaaaaaa|aaaaaaaa|aaabbbbb|bbbbbbbb|bbbbbb--|
     */
    DEBUG_ASSERT_LT(idataword, m_dsize_spin_channel, "dataword index OOB");

    const auto nbit_offset= m_nbit_in_last_alpha_dataword;
    const auto nbit_next = defs::nbit_word - nbit_offset;
    auto this_dword_mask = bit_utils::make_range_mask<size_t>(nbit_offset, defs::nbit_word);
    auto next_dword_mask = bit_utils::make_range_mask<size_t>(0, nbit_offset);
    std::cout << bit_utils::to_string(this_dword_mask) << std::endl;
    std::cout << bit_utils::to_string(next_dword_mask) << std::endl;
    std::cout << bit_utils::to_string(this_dword_mask >> nbit_offset) << std::endl;
    std::cout << bit_utils::to_string(next_dword_mask << nbit_next) << std::endl;

    auto dptr = reinterpret_cast<size_t *>(begin())+(m_dsize_spin_channel - 1);
    auto this_dword = (this_dword_mask & dptr[idataword]) >> nbit_offset;
    if (idataword+1==m_dsize_spin_channel){
        // last dword, there is no next_dword
        return this_dword;
    }
    auto next_dword = (next_dword_mask & dptr[idataword+1]) << nbit_next;
    return this_dword | next_dword;
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