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

void Determinant::set(const size_t &ispat, const size_t &ispin){
    assert(ispin==0 || ispin==1);
    assert(ispat>=0 && ispat<m_bitfields[0].m_nbit);
    m_bitfields[ispin].set(ispat);
}

void Determinant::set(const size_t &i){
    set(i%m_bitfields[0].m_nbit, i/m_bitfields[0].m_nbit);
}

void Determinant::set(const defs::inds &inds){
    for (auto ind: inds) set(ind);
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
    DeterminantSetEnumerator set_enumerator(*this);
    DeterminantXorEnumerator xor_enumerator(*this, other);

    bool counting_perms = false;
    size_t nperms = 0ul;

    size_t next_moved;
    xor_enumerator.next(next_moved);

    size_t iset;
    while(set_enumerator.next(iset)){
        if (iset>=next_moved){
            xor_enumerator.next(next_moved);
            counting_perms = !counting_perms;
        }
        else if (counting_perms) ++nperms;
    }
    return nperms%2;
}
