//
// Created by Robert J. Anderson on 08/08/2021.
//

#include "IntegralArray1e.h"

integrals_1e::IndexerSymNone::IndexerSymNone(size_t norb) : Indexer(norb, norb * norb) {}

size_t integrals_1e::IndexerSymNone::index_only(size_t a, size_t i) const {
    return a * m_norb + i;
}

std::pair<size_t, bool> integrals_1e::IndexerSymNone::index_and_conj(size_t a, size_t i) const {
    return {index_only(a, i), false};
}

void integrals_1e::IndexerSymNone::foreach(const std::function<void(size_t, size_t)> &fn) const {
    for (size_t a=0ul; a<m_norb; ++a){
        for (size_t i=0ul; i<m_norb; ++i) {
            fn(a, i);
        }
    }
}

integrals_1e::IndexerSymH::IndexerSymH(size_t norb) : Indexer(norb, npair(norb)) {}

size_t integrals_1e::IndexerSymH::index_only(size_t a, size_t i) const {
    return a>=i ? trigmap(a, i) : trigmap(i, a);
}

std::pair<size_t, bool> integrals_1e::IndexerSymH::index_and_conj(size_t a, size_t i) const {
    if (a>=i) return {trigmap(a, i), true};
    return {trigmap(i, a), false};
}

void integrals_1e::IndexerSymH::foreach(const std::function<void(size_t, size_t)> &fn) const {
    for (size_t a=0ul; a<m_norb; ++a){
        for (size_t i=a; i<m_norb; ++i) {
            fn(a, i);
        }
    }
}

std::string integrals_1e::syms::name(integrals_1e::syms::Sym sym) {
    switch (sym) {
        case Null: return "NULL";
        case None: return "no";
        case H: return "hermitian";
    }
    return {};
}
