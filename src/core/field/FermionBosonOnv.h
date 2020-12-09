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

        View(Field<FermionOnvSpecifier> &fonv_field, Field<BosonOnvSpecifier> &bonv_field,
             const size_t &irow, const size_t &iflat):
                m_fonv(fonv_field.get_view(irow, iflat)), m_bonv(bonv_field.get_view(irow, iflat)) {}


        View(const Field<FermionOnvSpecifier> &fonv_field, const Field<BosonOnvSpecifier> &bonv_field,
             const size_t &irow, const size_t &iflat):
                m_fonv(fonv_field.get_view(irow, iflat)), m_bonv(bonv_field.get_view(irow, iflat)) {}

        View &operator=(const std::pair<defs::inds, defs::inds> &pair) {
            m_fonv = pair.first;
            m_bonv = pair.second;
            return *this;
        }

        bool operator==(const View &other) const;

        bool operator!=(const View &other) const;

        void zero() {
            m_fonv.zero();
            m_bonv.zero();
        }

        bool is_zero() const {
            return m_fonv.is_zero() && m_bonv.is_zero();
        }

        void mpi_bcast(size_t iroot = 0) {
            m_fonv.mpi_bcast(iroot);
            m_bonv.mpi_bcast(iroot);
        }

        std::string to_string() const;

        void print() const;
    };


    template<size_t nind>
    struct Field : NdFieldGroup<nind> {

        NdFieldBase<FermionOnvSpecifier, nind> m_fonv;
        NdFieldBase<BosonOnvSpecifier, nind> m_bonv;

        using NdFieldGroup<nind>::m_format;

        template<typename ...Args>
        Field(Table *table, size_t nsite, std::string description, Args... shape) :
                NdFieldGroup<nind>(shape...),
                m_fonv(table, nsite, description + " (FermionOnv)", m_format),
                m_bonv(table, nsite, description + " (Boson ONV)", m_format) {}

        typedef View view_t;
        typedef const View const_view_t;

        template<typename ...Args>
        view_t operator()(const size_t &irow, Args... inds) {
            return view_t(m_fonv, m_bonv, irow, m_format.flatten(inds...));
        }

        template<typename ...Args>
        const_view_t operator()(const size_t &irow, Args... inds) const {
            return const_view_t(m_fonv, m_bonv, irow, m_format.flatten(inds...));
        }

        struct hash_fn {
            defs::hash_t operator()(const view_t &view) const {
                return FermionOnvSpecifier::hash(view.m_fonv) ^ BosonOnvSpecifier::hash(view.m_bonv);
            }
        };
    };

}


#endif //M7_FERMIONBOSONONV_H
