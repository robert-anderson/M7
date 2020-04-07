//
// Created by Robert John Anderson on 2020-03-30.
//

#include "Connection.h"
#include <algorithm>

Connection::Connection(const Field* field):
    m_element_dsize(field->element_dsize()), m_nbit(field->nbit()), m_ann(m_nbit), m_cre(m_nbit) {
}

Connection::Connection(const DeterminantElement &ket, const DeterminantElement &bra) : Connection(ket.field()) {
    assert(ket.compatible_with(bra));
    connect(ket, bra);
}

Connection::Connection(const DeterminantElement &ket) : Connection(ket, ket){}


void Connection::connect(const DeterminantElement &ket, const DeterminantElement &bra) {
    assert(ket.nbit() == m_nbit);
    assert(ket.dsize() == m_element_dsize);
    assert(ket.compatible_with(bra));
    m_nann = 0ul;
    m_ncre = 0ul;

    DeterminantElement::DatawordEnumerator ket_enumerator(ket);
    DeterminantElement::DatawordEnumerator bra_enumerator(bra);
    defs::data_t ket_work, bra_work, work;
    size_t idataword = ~0ul;
    while(ket_enumerator.next(ket_work, idataword) && bra_enumerator.next(bra_work)){
        work = ket_work &~ bra_work;
        while (work) m_ann[m_nann++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
        work = bra_work &~ ket_work;
        while (work) m_cre[m_ncre++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
    }
    assert(m_ncre<m_cre.size());
    assert(m_nann<m_ann.size());
}

void Connection::apply(const DeterminantElement &ket, DeterminantElement &bra){
    assert(m_ncre<m_cre.size());
    assert(m_nann<m_ann.size());
#ifndef NDEBUG
    for (size_t i=0ul; i<m_nann; ++i) assert(ket.get(m_ann[i]));
    for (size_t i=0ul; i<m_ncre; ++i) assert(!ket.get(m_cre[i]));
#endif
    bra = ket;
    for (size_t i=0ul; i<m_nann; ++i) bra.clr(m_ann[i]);
    for (size_t i=0ul; i<m_ncre; ++i) bra.set(m_cre[i]);
}

const size_t &Connection::nexcit() const {
    assert(m_ncre==m_nann);
    return m_ncre;
}


AntisymConnection::AntisymConnection(const Field *field): Connection(field), m_com(m_nbit){}

AntisymConnection::AntisymConnection(const DeterminantElement &ket, const DeterminantElement &bra) :
    Connection(ket, bra), m_com(m_nbit) {
    connect(ket, bra);
}

AntisymConnection::AntisymConnection(const DeterminantElement &ket) : AntisymConnection(ket, ket) {}

void AntisymConnection::connect(const DeterminantElement &ket, const DeterminantElement &bra) {
    Connection::connect(ket, bra);
    m_ncom = 0ul;
    size_t nperm = 0ul;

    DeterminantElement::DatawordEnumerator ket_enumerator(ket);
    DeterminantElement::DatawordEnumerator bra_enumerator(bra);
    defs::data_t ket_work, bra_work, work;
    size_t idataword = ~0ul;

    auto des_iter = m_ann.begin();
    const auto des_end = m_ann.begin() + m_nann;
    auto cre_iter = m_cre.begin();
    const auto cre_end = m_cre.begin() + m_ncre;

    while (ket_enumerator.next(ket_work, idataword) && bra_enumerator.next(bra_work)) {
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

void AntisymConnection::apply(const DeterminantElement &ket) {
    assert(m_ncre < m_nbit);
    assert(m_nann < m_nbit);
    assert(std::is_sorted(m_cre.begin(), m_cre.begin()+m_ncre));
    assert(std::is_sorted(m_ann.begin(), m_ann.begin()+m_nann));
    m_ncom = 0ul;
    size_t nperm = 0ul;

    DeterminantElement::DatawordEnumerator enumerator(ket);
    defs::data_t work;
    size_t idataword = ~0ul;

    auto des_iter = m_ann.begin();
    const auto des_end = m_ann.begin() + m_nann;
    auto cre_iter = m_cre.begin();
    const auto cre_end = m_cre.begin() + m_ncre;

    while (enumerator.next(work, idataword)){
        while (work) {
            auto &com = m_com[m_ncom];
            auto setbit = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
            if (des_iter != des_end && setbit==*des_iter) {
                des_iter++;
                nperm+=m_ncom;
                continue;
            }
            assert(cre_iter==cre_end || setbit!=*cre_iter);
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

void AntisymConnection::apply(const DeterminantElement &ket, DeterminantElement &bra) {
    apply(ket);
    Connection::apply(ket, bra);
}