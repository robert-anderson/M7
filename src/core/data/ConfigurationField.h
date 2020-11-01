//
// Created by RJA on 31/10/2020.
//

#ifndef M7_CONFIGURATIONFIELD_H
#define M7_CONFIGURATIONFIELD_H

#include "NdCompositeField.h"
#include "Fields.h"

template<size_t nind>
struct FermionBosonField : NdCompositeField<nind> {

    fields::Determinants<nind> det;
    fields::BosonOnvs<nind> perm;

    template<typename ...Args>
    FermionBosonField(TableX *table, size_t nsite, size_t nmode, NdFormat<nind> format) :
    NdCompositeField<nind>(table, format),
    det(this, {nsite}, "Determinant"),
    perm(this, {nmode}, "Boson ONV"){}

    struct View {
        DeterminantFieldX::view_t det;
        BosonOnvField::view_t perm;
    };

    typedef View view_t;
    typedef const View const_view_t;

    template<typename ...Args>
    view_t operator()(const size_t &irow, Args... inds) {
        return {det(irow, inds...), perm(irow, inds...)};
    }

    template<typename ...Args>
    const_view_t operator()(const size_t &irow, Args... inds) const {
        return {det(irow, inds...), perm(irow, inds...)};
    }

};


#endif //M7_CONFIGURATIONFIELD_H
