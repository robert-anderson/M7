//
// Created by rja on 08/08/2021.
//

#include "FrmOnvField.h"


FrmOnvField::FrmOnvField(Row *row, size_t nsite, std::string name) :
        base_t(row, {{2, nsite}, {"spin channel", "site"}}, name),
        m_nsite(nsite), m_nspinorb(m_format.m_nelement){
}

FrmOnvField::FrmOnvField(Row *row, BasisDims bd, std::string name) :
        FrmOnvField(row, bd.m_nsite, name){
    bd.require_pure_frm();
}

FrmOnvField::FrmOnvField(const FrmOnvField &other) :
        base_t(other), m_nsite(other.m_nsite), m_nspinorb(other.m_nspinorb){
}

FrmOnvField &FrmOnvField::operator=(std::pair<const defs::inds &, const defs::inds &> setbits) {
    // prezero the element
    zero();
    for (const auto &ind: setbits.first) set(ind);
    for (const auto &ind: setbits.second) set(ind+m_nsite);
    return *this;
}

void FrmOnvField::foreach(const std::function<void(const size_t &)> &body_fn) const {
    size_t work;
    for (size_t idataword = 0; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            body_fn(ibit);
        }
    }
}

void FrmOnvField::foreach_pair_inner(const size_t &ibit,
                                     const std::function<void(const size_t &, const size_t &)> &body_fn) const {
    size_t work;
    for (size_t idataword = 0; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t jbit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            if (jbit==ibit) return;
            body_fn(jbit, ibit);
        }
    }
}

void FrmOnvField::foreach_pair(const std::function<void(const size_t &, const size_t &)> &body_fn) const {
    size_t work;
    for (size_t idataword = 0; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            foreach_pair_inner(ibit, body_fn);
        }
    }
}

void FrmOnvField::foreach_pair(const std::function<void(const size_t &)> &body_fn_outer,
                               const std::function<void(const size_t &, const size_t &)> &body_fn_inner) const {
    size_t work;
    for (size_t idataword = 0; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            body_fn_outer(ibit);
            foreach_pair_inner(ibit, body_fn_inner);
        }
    }
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

size_t FrmOnvField::isite(const size_t &ibit) const {
    return isite(ibit, m_nsite);
}

size_t FrmOnvField::ispin(const size_t &ibit) const {
    return ispin(ibit, m_nsite);
}

size_t FrmOnvField::ibit(const size_t &ispin, const size_t &isite) const {
    return ibit(ispin, isite, m_nsite);
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
