//
// Created by RJA on 31/10/2020.
//

#ifndef M7_CONFIGURATIONFIELD_H
#define M7_CONFIGURATIONFIELD_H

#include "NdCompositeField.h"
#include "NdField.h"
#include "DeterminantField.h"
#include "BosonOnvField.h"

template<size_t nind>
struct FermionField : NdCompositeField<nind> {

    NdFieldX<DeterminantFieldX, nind> m_det;

    template<typename ...Args>
    FermionField(TableX *table, size_t nsite, std::string description, NdFormat<nind> format) :
            NdCompositeField<nind>(table, format),
            m_det(this, DeterminantFieldX(nsite), description + " (Determinant)") {}

    struct View : CompositeField::View {
        DeterminantFieldX::view_t m_det;

        View(DeterminantFieldX::view_t &&det) :
                m_det(std::move(det)) {}

        std::string to_string() const override {
            return m_det.to_string();
        }
    };

    typedef DeterminantFieldX::view_t view_t;
    typedef DeterminantFieldX::const_view_t const_view_t;

    template<typename ...Args>
    view_t operator()(const size_t &irow, Args... inds) {
        return m_det(irow, inds...);
    }

    template<typename ...Args>
    const_view_t operator()(const size_t &irow, Args... inds) const {
        return m_det(irow, inds...);
    }

};

template<size_t nind>
struct FermionBosonField : FermionField<nind> {

    using FermionField<nind>::m_det;
    NdFieldX<BosonOnvField, nind> m_perm;

    template<typename ...Args>
    FermionBosonField(TableX *table, size_t nsite, size_t nmode, std::string description, NdFormat<nind> format) :
            FermionField<nind>(table, nsite, description, format),
            m_perm(this, BosonOnvField(nmode), description + " (Boson ONV)") {}

    struct View : CompositeField::View {
        DeterminantFieldX::view_t m_det;
        BosonOnvField::view_t m_perm;

        View(DeterminantFieldX::view_t &&det,
             BosonOnvField::view_t &&perm) :
                m_det(std::move(det)), m_perm(std::move(perm)) {}

        bool operator==(const View &other) const {
            return m_det==other.m_det && m_perm==other.m_perm;
        }

        bool operator!=(const View &other) const {
            return !(*this==other);
        }

        std::string to_string() const override {
            return m_det.to_string() + " " + m_perm.to_string();
        }
    };

    typedef View view_t;
    typedef const View const_view_t;

    template<typename ...Args>
    view_t operator()(const size_t &irow, Args... inds) {
        return {m_det(irow, inds...), m_perm(irow, inds...)};
    }

    template<typename ...Args>
    const_view_t operator()(const size_t &irow, Args... inds) const {
        return {m_det(irow, inds...), m_perm(irow, inds...)};
    }

};


#endif //M7_CONFIGURATIONFIELD_H
