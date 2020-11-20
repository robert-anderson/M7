//
// Created by RJA on 31/10/2020.
//

#ifndef M7_FERMIONBOSONONV_H
#define M7_FERMIONBOSONONV_H

#include "FermionOnvSpecifier.h"
#include "BosonOnvSpecifier.h"
#include "TableField.h"

/*
 * Fermion-boson product state
 */

namespace fb_onv {
    struct View {
        FermionOnvSpecifier::view_t m_fonv;
        BosonOnvSpecifier::view_t m_bonv;

        View(FermionOnvSpecifier::view_t &&fonv,
             BosonOnvSpecifier::view_t &&bonv);

        View& operator=(const std::pair<defs::inds, defs::inds>& pair) {
            m_fonv = pair.first;
            m_bonv = pair.second;
            return *this;
        }

        bool operator==(const View &other) const;

        bool operator!=(const View &other) const;

        void zero() {
            m_fonv.zero(); m_bonv.zero();
        }

        bool is_zero() const {
            return m_fonv.is_zero() && m_bonv.is_zero();
        }

        std::string to_string();

        void print();
    };


    template<size_t nind>
    struct Field : NdFieldGroup<nind> {

        NdFieldBase<FermionOnvSpecifier, nind> m_fonv;
        NdFieldBase<BosonOnvSpecifier, nind> m_bonv;

        using NdFieldGroup<nind>::m_format;

        template<typename ...Args>
        Field(TableX *table, size_t nsite, std::string description, Args... shape) :
                NdFieldGroup<nind>(shape...),
                m_fonv(table, nsite, description + " (FermionOnv)", m_format),
                m_bonv(table, nsite, description + " (Boson ONV)", m_format) {}

        typedef View view_t;
        typedef const View const_view_t;

        template<typename ...Args>
        view_t operator()(const size_t &irow, Args... inds) {
            return {m_fonv(irow, inds...), m_bonv(irow, inds...)};
        }

        template<typename ...Args>
        const_view_t operator()(const size_t &irow, Args... inds) const {
            return {m_fonv(irow, inds...), m_bonv(irow, inds...)};
        }

        struct hash_fn {
            defs::hash_t operator()(const view_t &view) const {
                return FermionOnvSpecifier::hash(view.m_fonv) ^ BosonOnvSpecifier::hash(view.m_bonv);
            }
        };
    };

}


#endif //M7_FERMIONBOSONONV_H
