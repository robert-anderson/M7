//
// Created by Robert John Anderson on 2020-03-30.
//

#include "Connection.h"
#include <algorithm>

Connection::Connection(const DeterminantSpecifier& field):
    m_element_dsize(field.m_ndataword), m_nbit(field.m_nbit) {}

Connection::Connection(const views::Determinant &ket, const views::Determinant &bra) : Connection(ket.field()) {
    ASSERT(ket.nsite() == bra.nsite());
    connect(ket, bra);
}

Connection::Connection(const views::Determinant &ket) : Connection(ket, ket){}


void Connection::connect(const views::Determinant &ket, const views::Determinant &bra) {
    ASSERT(!ket.is_zero());
    ASSERT(!bra.is_zero());
    ASSERT(ket.nbit() == m_nbit);
    ASSERT(ket.ndataword() == m_element_dsize);
    ASSERT(ket.ndataword() == bra.ndataword());
    m_nann = 0ul;
    m_ncre = 0ul;


    defs::data_t ket_work, bra_work, work;
    for (size_t idataword = 0ul; idataword<ket.ndataword(); ++idataword){
        ket_work = ket.get_dataword(idataword);
        bra_work = ket.get_dataword(idataword);
        work = ket_work&~bra_work;
        while (work) m_ann[m_nann++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
        work = bra_work &~ ket_work;
        while (work) m_cre[m_ncre++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
    }
    ASSERT(m_ncre < m_cre.size());
    ASSERT(m_nann < m_ann.size());
}

void Connection::apply(const views::Determinant &ket, views::Determinant &bra){
    ASSERT(!ket.is_zero());
    ASSERT(m_ncre < m_cre.size());
    ASSERT(m_nann < m_ann.size());
#ifndef NDEBUG
    for (size_t i=0ul; i<m_nann; ++i) ASSERT(ket.get(m_ann[i]));
    for (size_t i=0ul; i<m_ncre; ++i) ASSERT(!ket.get(m_cre[i]));
#endif
    bra = ket;
    for (size_t i=0ul; i<m_nann; ++i) bra.clr(m_ann[i]);
    for (size_t i=0ul; i<m_ncre; ++i) bra.set(m_cre[i]);
    ASSERT(ket.nsetbit()==bra.nsetbit());
}

const size_t &Connection::nexcit() const {
    ASSERT(m_ncre == m_nann);
    return m_ncre;
}


AntisymConnection::AntisymConnection(const DeterminantSpecifier &field): Connection(field) {}

AntisymConnection::AntisymConnection(const views::Determinant &ket, const views::Determinant &bra) :
    Connection(ket, bra) {
    connect(ket, bra);
}

AntisymConnection::AntisymConnection(const views::Determinant &ket) : AntisymConnection(ket.field()) {}

void AntisymConnection::connect(const views::Determinant &ket, const views::Determinant &bra) {
    Connection::connect(ket, bra);
    m_ncom = 0ul;
    size_t nperm = 0ul;

    auto des_iter = m_ann.begin();
    const auto des_end = m_ann.begin() + m_nann;
    auto cre_iter = m_cre.begin();
    const auto cre_end = m_cre.begin() + m_ncre;

    defs::data_t ket_work, bra_work, work;
    for (size_t idataword = 0ul; idataword<ket.ndataword(); ++idataword){
        ket_work = ket.get_dataword(idataword);
        bra_work = ket.get_dataword(idataword);
        work = ket_work & bra_work;
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

void AntisymConnection::apply(const views::Determinant &ket) {
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

    defs::data_t work;
    for(size_t idataword=0ul; idataword<ket.ndataword(); ++idataword){
        work = ket.get_dataword(idataword);
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

void AntisymConnection::apply(const views::Determinant &ket, views::Determinant &bra) {
    apply(ket);
    Connection::apply(ket, bra);
}
