//
// Created by Robert John Anderson on 2020-03-30.
//

#include "FermionOnvConnection.h"
#include <algorithm>

FermionOnvConnection::FermionOnvConnection(const FermionOnvSpecifier& spec):
    m_element_dsize(spec.m_ndataword), m_nbit(spec.m_nbit) {}

FermionOnvConnection::FermionOnvConnection(const views::FermionOnv &in, const views::FermionOnv &out) : FermionOnvConnection(in.spec()) {
    ASSERT(in.nsite() == out.nsite());
    connect(in, out);
}

FermionOnvConnection::FermionOnvConnection(const views::FermionOnv &in) : FermionOnvConnection(in, in){}


void FermionOnvConnection::connect(const views::FermionOnv &in, const views::FermionOnv &out) {
    ASSERT(!in.is_zero());
    ASSERT(!out.is_zero());
    ASSERT(in.nbit() == m_nbit);
    ASSERT(in.ndataword() == m_element_dsize);
    ASSERT(in.ndataword() == out.ndataword());
    m_nann = 0ul;
    m_ncre = 0ul;


    defs::data_t in_work, out_work, work;
    for (size_t idataword = 0ul; idataword<in.ndataword(); ++idataword){
        in_work = in.get_dataword(idataword);
        out_work = out.get_dataword(idataword);
        work = in_work&~out_work;
        while (work) m_ann[m_nann++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
        work = out_work &~ in_work;
        while (work) m_cre[m_ncre++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
    }
    ASSERT(m_ncre < m_cre.size());
    ASSERT(m_nann < m_ann.size());
}

void FermionOnvConnection::apply(const views::FermionOnv &in, views::FermionOnv &out){
    ASSERT(!in.is_zero());
    ASSERT(m_ncre < m_cre.size());
    ASSERT(m_nann < m_ann.size());
#ifndef NDEBUG
    for (size_t i=0ul; i<m_nann; ++i) ASSERT(in.get(m_ann[i]));
    for (size_t i=0ul; i<m_ncre; ++i) ASSERT(!in.get(m_cre[i]));
#endif
    out = in;
    for (size_t i=0ul; i<m_nann; ++i) out.clr(m_ann[i]);
    for (size_t i=0ul; i<m_ncre; ++i) out.set(m_cre[i]);
    ASSERT(in.nsetbit()==out.nsetbit());
}

const size_t &FermionOnvConnection::nexcit() const {
    ASSERT(m_ncre == m_nann);
    return m_ncre;
}


AntisymFermionOnvConnection::AntisymFermionOnvConnection(const FermionOnvSpecifier &field): FermionOnvConnection(field) {}

AntisymFermionOnvConnection::AntisymFermionOnvConnection(const views::FermionOnv &in, const views::FermionOnv &out) :
        FermionOnvConnection(in, out) {
    connect(in, out);
}

AntisymFermionOnvConnection::AntisymFermionOnvConnection(const views::FermionOnv &in) : AntisymFermionOnvConnection(in.spec()) {}

void AntisymFermionOnvConnection::connect(const views::FermionOnv &in, const views::FermionOnv &out) {
    FermionOnvConnection::connect(in, out);
    m_ncom = 0ul;
    size_t nperm = 0ul;

    auto des_iter = m_ann.begin();
    const auto des_end = m_ann.begin() + m_nann;
    auto cre_iter = m_cre.begin();
    const auto cre_end = m_cre.begin() + m_ncre;

    defs::data_t in_work, out_work, work;
    for (size_t idataword = 0ul; idataword<in.ndataword(); ++idataword){
        in_work = in.get_dataword(idataword);
        out_work = out.get_dataword(idataword);
        work = in_work & out_work;
        while (work) {
            auto &com = m_com[m_ncom];
            com = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
            while (des_iter != des_end && *des_iter < com) {
                // an annihilation operator has been passed in the iteration over common indices
                des_iter++;
                nperm += m_ncom;
            }
            while (cre_iter != cre_end && *cre_iter < com) {
                // a creation operator has been passed in the iteration over common indices
                cre_iter++;
                nperm += m_ncom;
            }
            m_ncom++;
        }
    }
    while (des_iter != des_end) {des_iter++; nperm += m_ncom;}
    while (cre_iter != cre_end) {cre_iter++; nperm += m_ncom;}
    m_phase = nperm & 1ul;
}

void AntisymFermionOnvConnection::apply(const views::FermionOnv &in) {
    ASSERT(m_ncre < m_nbit);
    ASSERT(m_nann < m_nbit);
    ASSERT(std::is_sorted(m_cre.begin(), m_cre.begin() + m_ncre));
    ASSERT(std::is_sorted(m_ann.begin(), m_ann.begin() + m_nann));
    m_ncom = 0ul;
    size_t nperm = 0ul;

    auto des_iter = m_ann.begin();
    const auto des_end = m_ann.begin() + m_nann;
    auto cre_iter = m_cre.begin();
    const auto cre_end = m_cre.begin() + m_ncre;

    for(size_t idataword=0ul; idataword<in.ndataword(); ++idataword){
        auto work = in.get_dataword(idataword);
        while (work) {
            auto &com = m_com[m_ncom];
            auto setbit = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
            if (des_iter != des_end && setbit==*des_iter) {
                des_iter++;
                nperm+=m_ncom;
                continue;
            }
            // check we aren't trying to create an electron in an occupied orbital
            ASSERT((cre_iter == cre_end) || (setbit != *cre_iter));
            com = setbit;
            while (cre_iter != cre_end && *cre_iter < com) {
                cre_iter++;
                nperm += m_ncom;
            }
            m_ncom++;
        }
    }
    while (cre_iter != cre_end) {cre_iter++; nperm += m_ncom;}
    m_phase = nperm & 1ul;
}

void AntisymFermionOnvConnection::apply(const views::FermionOnv &in, views::FermionOnv &out) {
    apply(in);
    FermionOnvConnection::apply(in, out);
}

void AntisymFermionOnvConnection::zero() {
    FermionOnvConnection::zero();
    m_ncom = 0;
}
