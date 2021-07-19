//
// Created by Robert John Anderson on 2020-03-30.
//

#include "FermionOnvConnection.h"
#include <algorithm>

FermionOnvConnection::FermionOnvConnection(size_t nsite){
    m_ann.reserve(2*nsite);
    m_cre.reserve(2*nsite);
}

FermionOnvConnection::FermionOnvConnection(const fields::Onv<0> &in, const fields::Onv<0> &out):
        FermionOnvConnection(in.m_nsite) {
    ASSERT(in.m_nsite == out.m_nsite);
    connect(in, out);
}

FermionOnvConnection::FermionOnvConnection(const fields::Onv<0> &in):
        FermionOnvConnection(in, in){}

void FermionOnvConnection::connect(const fields::Onv<0> &in, const fields::Onv<0> &out) {
    ASSERT(!in.is_zero())
    zero();

    size_t in_work, out_work, work;
    for (size_t idataword = 0ul; idataword<in.m_dsize; ++idataword){
        in_work = in.get_dataword(idataword);
        out_work = out.get_dataword(idataword);
        work = in_work&~out_work;
        while (work) add_ann(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
        work = out_work &~ in_work;
        while (work) add_cre(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
    }
}

void FermionOnvConnection::apply(const fields::Onv<0> &in, fields::Onv<0> &out){
    ASSERT(!in.is_zero());
#ifndef NDEBUG
    for (size_t i=0ul; i<nann(); ++i) ASSERT(in.get(ann(i)));
    for (size_t i=0ul; i<ncre(); ++i) ASSERT(!in.get(cre(i)));
#endif
    out = in;
    for (size_t i=0ul; i<nann(); ++i) out.clr(ann(i));
    for (size_t i=0ul; i<ncre(); ++i) out.set(cre(i));
    ASSERT(in.nsetbit()==out.nsetbit());
}

size_t FermionOnvConnection::nexcit() const {
    ASSERT(nann() == ncre());
    return ncre();
}


AntisymFermionOnvConnection::AntisymFermionOnvConnection(size_t nsite):
FermionOnvConnection(nsite) {
    m_com.reserve(2*nsite);
}

AntisymFermionOnvConnection::AntisymFermionOnvConnection(const fields::Onv<0> &in, const fields::Onv<0> &out) :
        AntisymFermionOnvConnection(in) {
    connect(in, out);
}

AntisymFermionOnvConnection::AntisymFermionOnvConnection(const fields::Onv<0> &in) : AntisymFermionOnvConnection(in.m_nsite) {}

void AntisymFermionOnvConnection::connect(const fields::Onv<0> &in, const fields::Onv<0> &out) {
    FermionOnvConnection::connect(in, out);
    m_com.clear();
    size_t nperm = 0ul;

    auto ann_iter = m_ann.begin();
    auto cre_iter = m_cre.begin();

    size_t in_work, out_work, work;
    for (size_t idataword = 0ul; idataword<in.m_dsize; ++idataword){
        in_work = in.get_dataword(idataword);
        out_work = out.get_dataword(idataword);
        work = in_work & out_work;
        while (work) {
            auto setbit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            while (ann_iter != m_ann.end() && *ann_iter < setbit) {
                // an annihilation operator has been passed in the iteration over common indices
                ann_iter++;
                nperm += ncom();
            }
            while (cre_iter != m_cre.end() && *cre_iter < setbit) {
                // a creation operator has been passed in the iteration over common indices
                cre_iter++;
                nperm += ncom();
            }
            m_com.push_back(setbit);
        }
    }
    while (ann_iter != m_ann.end()) {ann_iter++; nperm += ncom();}
    while (cre_iter != m_cre.end()) {cre_iter++; nperm += ncom();}
    m_phase = nperm & 1ul;
}

void AntisymFermionOnvConnection::apply(const fields::Onv<0> &in) {
    ASSERT(std::is_sorted(m_cre.begin(), m_cre.end()));
    ASSERT(std::is_sorted(m_ann.begin(), m_ann.end()));
    m_com.clear();
    size_t nperm = 0ul;

    auto ann_iter = m_ann.begin();
    auto cre_iter = m_cre.begin();

    for(size_t idataword=0ul; idataword<in.m_dsize; ++idataword){
        auto work = in.get_dataword(idataword);
        while (work) {
            auto setbit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            if (ann_iter != m_ann.end() && setbit==*ann_iter) {
                ann_iter++;
                nperm+=ncom();
                continue;
            }
            // check we aren't trying to create an electron in an occupied orbital
            ASSERT((cre_iter == m_cre.end()) || (setbit != *cre_iter));
            while (cre_iter != m_cre.end() && *cre_iter < setbit) {
                cre_iter++;
                nperm += ncom();
            }
            m_com.push_back(setbit);
        }
    }
    while (cre_iter != m_cre.end()) {
        cre_iter++;
        nperm += ncom();
    }
    m_phase = nperm & 1ul;
}

void AntisymFermionOnvConnection::apply(const fields::Onv<0> &in, fields::Onv<0> &out) {
    apply(in);
    FermionOnvConnection::apply(in, out);
}

void AntisymFermionOnvConnection::zero() {
    FermionOnvConnection::zero();
    m_com.clear();
}