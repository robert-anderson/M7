//
// Created by Robert J. Anderson on 1/30/22.
//

#ifndef M7_XONVFIELD_H
#define M7_XONVFIELD_H

#include <M7_lib/basis/BasisData.h>

#include "CompositeField.h"
#include "Row.h"
#include "FrmOnvField.h"
#include "BosOnvField.h"
#include "FrmBosOnvField.h"

template<typename T>
struct XonvField : CompositeField<T, T> {
    using CompositeFieldBase::prefix;
    T m_ket, m_bra;
    /*
     * all ONVs implement a ctor parametrised by sys::Basis
     */
    XonvField(Row *row, const sys::Basis& basis, std::string name) :
            CompositeField<T, T>(m_ket, m_bra),
            m_ket(row, basis, prefix("ket", name)),
            m_bra(row, basis, prefix("bra", name)) {}

    XonvField(Row *row, const sys::Sector& sector, std::string name) : XonvField(row, sector.basis(), name){}
};


struct FrmXonvField : XonvField<FrmOnvField> {
    FrmXonvField(Row *row, const sys::Basis& basis, std::string name = ""): XonvField<FrmOnvField>(row, basis, name){}
    FrmXonvField(Row *row, const sys::frm::Basis& basis, std::string name = ""):
            FrmXonvField(row, sys::Basis(basis, {0ul}), std::move(name)){}
    FrmXonvField(Row *row, const sys::Sector& sector, std::string name = "") :
            FrmXonvField(row, sector.basis(), std::move(name)){}

    FrmXonvField &operator=(const std::pair<uintv_t, uintv_t> &inds) {
        m_ket = inds.first;
        m_bra = inds.second;
        return *this;
    }
};


struct BosXonvField : XonvField<BosOnvField> {

    BosXonvField(Row *row, const sys::Basis& basis, std::string name = ""): XonvField<BosOnvField>(row, basis, name){}
    BosXonvField(Row *row, const sys::bos::Basis& basis, std::string name = ""):
            BosXonvField(row, sys::Basis({0ul}, basis), std::move(name)){}
    BosXonvField(Row *row, const sys::Sector& sector, std::string name = "") :
            BosXonvField(row, sector.basis(), std::move(name)){}

    BosXonvField &operator=(const std::pair<uintv_t, uintv_t> &inds) {
        m_ket = inds.first;
        m_bra = inds.second;
        return *this;
    }
};

struct FrmBosXonvField : XonvField<FrmBosOnvField> {
    FrmBosXonvField(Row *row, const sys::Basis& basis, std::string name = ""):
            XonvField<FrmBosOnvField>(row, basis, std::move(name)){}

    FrmBosXonvField(Row *row, const sys::frm::Basis& frm_basis, const sys::bos::Basis& bos_basis, std::string name = ""):
            XonvField<FrmBosOnvField>(row, sys::Basis(frm_basis, bos_basis), std::move(name)){}

    FrmBosXonvField(Row *row, const sys::Sector& sector, std::string name = "") :
            FrmBosXonvField(row, sector.basis(), std::move(name)){}

    typedef std::pair<uintv_t, uintv_t> pair_t;
    FrmBosXonvField& operator=(const std::pair<pair_t, pair_t> &inds) {
        m_ket = inds.first;
        m_bra = inds.second;
        return *this;
    }
};

#endif //M7_XONVFIELD_H
