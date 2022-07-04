//
// Created by Robert J. Anderson on 08/08/2021.
//

#include "IntegralArray2e.h"
#include "M7_lib/foreach/BasicForeach.h"


str_t integrals_2e::syms::name(integrals_2e::syms::Sym sym) {
    switch (sym) {
        case Null: return "NULL";
        case None: return "no";
        case H: return "2-fold (hermiticity)";
        case D: return "2-fold (dummy integration variable interchange)";
        case DH: return "4-fold (complex orbitals)";
        case DR: return "4-fold (non-hermitian)";
        case DHR: return "8-fold";
    }
    return {};
}

strv_t integrals_2e::syms::equivalences(integrals_2e::syms::Sym sym) {
    switch (sym) {
        case Null: return {};
        case None: return {"<ab|ij>"};
        case H: return {"<ab|ij>", "<ij|ab>"};
        case D: return {"<ab|ij>", "<ba|ji>"};
        case DH: return {"<ab|ij>", "<ba|ji>", "<ij|ab>*", "<ji|ba>*"};
        case DR: return {"<ab|ij>", "<aj|ib>", "<ba|ji>", "<ja|bi>"};
        case DHR: return {"<ab|ij>", "<aj|ib>", "<ba|ji>", "<ja|bi>",
                          "<aj|ib>", "<ab|ij>", "<bi|ja>", "<ji|ba>"};
    }
    return {};
}

integrals_2e::IndexerSymNone::IndexerSymNone(uint_t norb) :
        Indexer(norb, pow(norb, 4ul), syms::None), m_norb2(norb * norb), m_norb3(norb * m_norb2) {}

uint_t integrals_2e::IndexerSymNone::index_only(uint_t a, uint_t b, uint_t i, uint_t j) const {
    return a * m_norb3 + b * m_norb2 + i * m_norb + j;
}

std::pair<uint_t, bool> integrals_2e::IndexerSymNone::index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const {
    return {index_only(a, b, i, j), false};
}

void integrals_2e::Indexer::foreach(const integrals_2e::foreach_fn_t &fn) const {
    using namespace basic_foreach::ctnd;
    Unrestricted<4> foreach(m_norb);
    foreach.loop([&fn](const inds_t<4>& inds){fn(inds[0], inds[1], inds[2], inds[3]);});
}


integrals_2e::IndexerSymH::IndexerSymH(uint_t norb) : Indexer(norb, npair(norb*norb), syms::H) {}

uint_t integrals_2e::IndexerSymH::index_only(uint_t a, uint_t b, uint_t i, uint_t j) const {
    const auto ab = a*m_norb+b;
    const auto ij = i*m_norb+j;
    return (ab>=ij) ? trigmap(ab, ij) : trigmap(ij, ab);
}

std::pair<uint_t, bool> integrals_2e::IndexerSymH::index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const {
    const auto ab = a*m_norb+b;
    const auto ij = i*m_norb+j;
    if (ab>=ij) return {trigmap(ab, ij), true};
    else return {trigmap(ij, ab), false};
}

uint_t integrals_2e::IndexerSymD::index_only(uint_t a, uint_t b, uint_t i, uint_t j) const {
    return m_sym_h.index_only(a, i, b, j);
}

integrals_2e::IndexerSymDH::IndexerSymDH(uint_t norb) :
        Indexer(norb, 2*npair(npair(norb)), syms::DH), m_hir_size(m_size/2){}

std::pair<uint_t, bool> integrals_2e::IndexerSymD::index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const {
    return {index_only(a, b, i, j), false};
}

integrals_2e::IndexerSymD::IndexerSymD(uint_t norb) : Indexer(norb, npair(norb*norb), syms::D), m_sym_h(norb){}

std::pair<uint_t, bool> integrals_2e::IndexerSymDH::index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const {
    // let y = <ab|ij>, correct order is x
    if (b >= j){
        const auto bj = trigmap(b, j);
        if (a >= i) {
            const auto ai= trigmap(a, i);
            if (bj >= ai) {
                // x = y
                return {trigmap(bj, ai), false};
            }
            else {
                // x = Dy i.e. y = Dx
                return {trigmap(ai, bj), false};
            }
        }
        else {
            const auto ai= trigmap(i, a);
            if (bj >= ai) {
                // x = KHRy i.e. y = RHKx = HRKx
                return {trigmap(bj, ai)+m_hir_size, true};
            }
            else {
                // x = DKHRy i.e. y = RHDKx = DRx
                return {trigmap(ai, bj)+m_hir_size, false};
            }
        }
    }
    else {
        const auto bj = trigmap(j, b);
        if (a >= i){
            const auto ai= trigmap(a, i);
            if (bj>=ai) {
                // x = Ry i.e. y = Rx
                return {trigmap(bj, ai)+m_hir_size, false};
            }
            else {
                // x = DRy i.e. y = RDx = HDRKx
                return {trigmap(ai, bj)+m_hir_size, true};
            }
        }
        else {
            const auto ai= trigmap(i, a);
            if (bj>=ai) {
                // x = KHRRy i.e. y = HKx
                return {trigmap(bj, ai), true};
            }
            else {
                // x = DKHRRy i.e. y = HDKx
                return {trigmap(ai, bj), true};
            }
        }
    }
    return {};
}

uint_t integrals_2e::IndexerSymDH::index_only(uint_t a, uint_t b, uint_t i, uint_t j) const {
    return index_and_conj(a, b, i, j).first;
}

integrals_2e::IndexerSymDR::IndexerSymDR(uint_t norb) :
        Indexer(norb, 2*npair(npair(norb)), syms::DR), m_hir_size(m_size/2){}

uint_t integrals_2e::IndexerSymDR::index_only(uint_t a, uint_t b, uint_t i, uint_t j) const {
    // let y = <ab|ij>, correct order is x
    if (b >= j){
        const auto bj = trigmap(b, j);
        if (a >= i) {
            const auto ai= trigmap(a, i);
            if (bj >= ai) {
                // x = y
                return trigmap(bj, ai);
            }
            else {
                // x = Dy i.e. y = Dx
                return trigmap(ai, bj);
            }
        }
        else {
            const auto ai= trigmap(i, a);
            if (bj >= ai) {
                // x = KHRy i.e. y = RHKx = HRKx
                return trigmap(bj, ai)+m_hir_size;
            }
            else {
                // x = DKHRy i.e. y = RHDKx = DRx
                return trigmap(ai, bj);
            }
        }
    }
    else {
        const auto bj = trigmap(j, b);
        if (a >= i){
            const auto ai= trigmap(a, i);
            if (bj>=ai) {
                // x = Ry i.e. y = Rx
                return trigmap(bj, ai);
            }
            else {
                // x = DRy i.e. y = RDx = HDRKx
                return trigmap(ai, bj)+m_hir_size;
            }
        }
        else {
            const auto ai= trigmap(i, a);
            if (bj>=ai) {
                // x = KHRRy i.e. y = HKx
                return trigmap(bj, ai)+m_hir_size;
            }
            else {
                // x = DKHRRy i.e. y = HDKx
                return trigmap(ai, bj)+m_hir_size;
            }
        }
    }
    return {};
}

std::pair<uint_t, bool> integrals_2e::IndexerSymDR::index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const {
    return {index_only(a, b, i, j), false};
}

integrals_2e::IndexerSymDHR::IndexerSymDHR(uint_t norb) : Indexer(norb, npair(npair(norb)), syms::DHR){}

uint_t integrals_2e::IndexerSymDHR::index_only(uint_t a, uint_t b, uint_t i, uint_t j) const {
    const auto bj = (b>=j) ? trigmap(b, j) : trigmap(j, b);
    const auto ai = (a>=i) ? trigmap(a, i) : trigmap(i, a);
    return (bj>=ai) ? trigmap(bj, ai) : trigmap(ai, bj);
}

std::pair<uint_t, bool> integrals_2e::IndexerSymDHR::index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const {
    // always have real orbitals so never any conjugation
    return {index_only(a, b, i, j), false};
}
