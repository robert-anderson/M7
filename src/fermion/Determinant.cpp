//
// Created by Robert John Anderson on 2020-01-27.
//

#include <iostream>
#include <assert.h>
#include "Determinant.h"
#include "../enumerators/BitfieldEnumerator.h"

Determinant::Determinant(const size_t &nspatorb) :
m_bitfields({BitfieldNew(nspatorb), BitfieldNew(nspatorb)}){}

Determinant::Determinant(const size_t &nspatorb, defs::data_t* data1, defs::data_t* data2) :
m_bitfields({BitfieldNew(nspatorb, data1), BitfieldNew(nspatorb, data2)}){}

Determinant::Determinant(const BitfieldNew &data1, const BitfieldNew &data2):
    m_bitfields({data1, data2}){}

std::string Determinant::to_string() const {
    return m_bitfields[0].to_string()+" "+m_bitfields[1].to_string();
}

void Determinant::print() const {
    std::cout << to_string() << std::endl;
}

size_t Determinant::nexcit(const Determinant &other) const {
    return m_bitfields[0].nsetbits_cleared(other.m_bitfields[0])+
           m_bitfields[1].nsetbits_cleared(other.m_bitfields[1]);
}

void Determinant::zero() {
    m_bitfields[0].zero();
    m_bitfields[1].zero();
}

bool Determinant::is_zero() const {
    return m_bitfields[0].is_zero() && m_bitfields[1].is_zero();
}


void Determinant::set(const size_t &ispat, const size_t &ispin){
    assert(ispin==0 || ispin==1);
    assert(ispat>=0 && ispat<m_bitfields[0].m_nbit);
    m_bitfields[ispin].set(ispat);
}

void Determinant::set(const size_t &i){
    assert(i<nspatorb()*2);
    set(i%m_bitfields[0].m_nbit, i/m_bitfields[0].m_nbit);
}

void Determinant::set(const defs::inds &inds){
    zero();
    for (auto ind: inds) set(ind);
}

void Determinant::clr(const size_t &ispat, const size_t &ispin){
    assert(ispin==0 || ispin==1);
    assert(ispat>=0 && ispat<m_bitfields[0].m_nbit);
    m_bitfields[ispin].clr(ispat);
}

void Determinant::clr(const size_t &i){
    assert(i<nspatorb()*2);
    clr(i%m_bitfields[0].m_nbit, i/m_bitfields[0].m_nbit);
}

bool Determinant::partial_phase(const defs::inds &removed, const size_t &nremoved) const{
    size_t nperm = 0ul;
    DeterminantSetEnumerator set_enumerator(*this);
    size_t set;
    for (auto iorb : removed){
        while (set_enumerator.next(set) && set < iorb) nperm++;
    }
    return nperm&1ul;
}

bool Determinant::partial_phase(const Determinant &other) const{
    size_t nperm = 0ul;
    size_t nperm_tot = 0ul;
    DeterminantSetEnumerator set_enumerator(*this);
    DeterminantAndNotEnumerator nand_enumerator(*this, other);
    size_t removed;
    size_t set;
    while(nand_enumerator.next(removed)){
        while (set_enumerator.next(set) && set < removed) nperm++;
        nperm_tot+=nperm;
    }
    return nperm_tot&1ul;
}

bool Determinant::phase(const Determinant &other) const {

    /*
     * Two determinants bra and ket are connected by ordered arrays of:
     *   -  nremoved spin-orbital indices occupied in the ket but not the bra
     *   -  ninserted spin-orbital indices occupied in the bra but not the ket
     * called "removed" and "inserted" respectively.
     * nremoved and ninserted need not be equal.
     *
     * let the removed and inserted spin orbitals be denoted R_i and I_i
     *
     *       __ ninserted-1  ^+    __ nremoved-1  ^
     * <bra| ||              a     ||             a   |ket> = +/- 1
     *          i=0           I_i     i=0          R_i
     */
    return partial_phase(other)^other.partial_phase(*this);
}

size_t Determinant::nelec() const {
    return m_bitfields[0].nsetbits()+m_bitfields[1].nsetbits();
}

size_t Determinant::nspatorb() const {
    return m_bitfields[0].m_nbit;
}

bool Determinant::operator==(const Determinant &rhs) const {
    return compare(rhs)==0;
}

int Determinant::compare(const Determinant &rhs) const {
    auto t0 = m_bitfields[0].compare(rhs.m_bitfields[0]);
    auto t1 = m_bitfields[1].compare(rhs.m_bitfields[1]);
    return t0?t0:t1;
}

Determinant Determinant::get_excited_det(const defs::inds &removed, const defs::inds &inserted) const{
    auto out = *this;
    for (auto iremoved : removed) out.clr(iremoved);
    for (auto iinserted : inserted) out.set(iinserted);
    return out;
}

Determinant Determinant::get_excited_det(const size_t &removed, const size_t &inserted) const{
    auto out = *this;
    out.clr(removed);
    out.set(inserted);
    return out;
}


Determinant Determinant::get_excited_det(
        const size_t &removed1, const size_t &removed2,
        const size_t &inserted1, const size_t &inserted2
) const{
    auto out = *this;
    out.clr(removed1);
    out.clr(removed2);
    out.set(inserted1);
    out.set(inserted2);
    return out;
}

Determinant &Determinant::operator=(const Determinant &rhs) {
    m_bitfields[0] = rhs.m_bitfields[0];
    m_bitfields[1] = rhs.m_bitfields[1];
    return *this;
}

bool Determinant::is_null() const {
    return nspatorb()==0;
}
