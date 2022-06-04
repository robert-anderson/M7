//
// Created by Robert J. Anderson on 08/08/2021.
//

#include "IntegralArray1e.h"

integrals_1e::IndexerSymNone::IndexerSymNone(size_t norb) : IntegralIndexer(norb, norb * norb) {}

size_t integrals_1e::IndexerSymNone::index_only(size_t a, size_t i) const {
    return a * m_norb + i;
}

std::pair<size_t, bool> integrals_1e::IndexerSymNone::index_and_conj(size_t a, size_t i) const {
    return {index_only(a, i), false};
}

integrals_1e::IndexerSymH::IndexerSymH(size_t norb) : IntegralIndexer(norb, npair(norb)) {}

size_t integrals_1e::IndexerSymH::index_only(size_t a, size_t i) const {
    return a>=i ? trigmap(a, i) : trigmap(i, a);
}

std::pair<size_t, bool> integrals_1e::IndexerSymH::index_and_conj(size_t a, size_t i) const {
    if (a>=i) return {trigmap(a, i), true};
    return {trigmap(i, a), false};
}
