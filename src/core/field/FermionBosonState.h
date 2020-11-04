//
// Created by RJA on 31/10/2020.
//

#ifndef M7_FERMIONBOSONSTATE_H
#define M7_FERMIONBOSONSTATE_H

#include "DeterminantSpecifier.h"
#include "BosonOnvSpecifier.h"
#include "TableField.h"

/*
 * Fermion-boson product state
 */

namespace fb_state {
    struct View {
        DeterminantSpecifier::view_t m_det;
        BosonOnvSpecifier::view_t m_perm;

        View(DeterminantSpecifier::view_t &&det,
             BosonOnvSpecifier::view_t &&perm) :
                m_det(std::move(det)), m_perm(std::move(perm)) {}

        bool operator==(const View &other) const {
            return m_det==other.m_det && m_perm==other.m_perm;
        }

        bool operator!=(const View &other) const {
            return !(*this==other);
        }

        std::string to_string() {
            return m_det.to_string() + " " + m_perm.to_string();
        }

        void print() {
            std::cout << to_string() << std::endl;
        }
    };


    template<size_t nind>
    struct Field : NdFieldGroup<nind> {

        NdFieldBase<DeterminantSpecifier, nind> m_det;
        NdFieldBase<BosonOnvSpecifier, nind> m_perm;

        using NdFieldGroup<nind>::m_format;
        template<typename ...Args>
        Field(TableX *table, size_t nsite, size_t nmode, std::string description, Args... shape) :
        NdFieldGroup<nind>(shape...),
        m_det(table, {nsite}, description + " (Determinant)", m_format),
        m_perm(table, {nmode}, description + " (Boson ONV)", m_format){}

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

}




#endif //M7_FERMIONBOSONSTATE_H
