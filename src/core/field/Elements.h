//
// Created by rja on 01/11/2020.
//

#ifndef M7_ELEMENTS_H
#define M7_ELEMENTS_H

#include "src/core/table/Table.h"
#include "src/core/table/BufferedTable.h"
#include "Fields.h"

template<typename viewable_t>
struct SingleFieldTable: TableX {
    static_assert(std::is_base_of<NdFieldGroup<0ul>, viewable_t>::value, "Template arg must be a scalar NdFieldGroup");
    viewable_t m_field;
    template<typename ...Args>
    SingleFieldTable(Args... args):
    TableX(),
    m_field(this, args...){}
};

template<typename viewable_t>
struct BufferedSingleFieldTable: BufferedTable<SingleFieldTable<viewable_t>>{
    typedef BufferedTable<SingleFieldTable<viewable_t>> base_t;
    template<typename ...Args>
    BufferedSingleFieldTable(Args... args): base_t(args...){
        base_t::expand(1ul);
        base_t::push_back();
    }
};

template <typename viewable_t>
struct Element : BufferedSingleFieldTable<viewable_t>, viewable_t::view_t {
    typedef BufferedSingleFieldTable<viewable_t> base_t;
    template<typename ...Args>
    Element(Args... args): base_t(args...), viewable_t::view_t(base_t::m_field(0)){}
};

namespace elements {
    struct Determinant : Element<fields::FermionOnv> {
        Determinant(size_t nsite):
        Element<fields::FermionOnv>(FermionOnvSpecifier(nsite), "Working determinant"){}
        using specs::FermionOnv::view_t::operator=;
    };

    struct BosonOnv : Element<fields::BosonOnv> {
        BosonOnv(size_t nmode):
                Element<fields::BosonOnv>(BosonOnvSpecifier(nmode), "Working boson ONV"){}
        using specs::BosonOnv::view_t::operator=;
    };

    struct FermionBosonConfiguration : Element<fields::FermiBosOnv> {
        FermionBosonConfiguration(size_t nsite, size_t nmode):
        Element<fields::FermiBosOnv>(nsite, nmode, "Working fermion-boson configuration"){}
    };

    template<size_t nind, bool bosons>
    struct ConfigurationSelector {};

    template<size_t nind>
    struct ConfigurationSelector<nind, false> {
        typedef Determinant type;
    };

    template<size_t nind>
    struct ConfigurationSelector<nind, true> {
        typedef FermionBosonConfiguration type;
    };

    using Configuration = ConfigurationSelector<0ul, defs::bosons>::type;
}
#if 0
#include "BufferedField.h"
#include "BufferedComposite.h"
#include "NumericField.h"
#include "NumericArraySpecifier.h"
#include "Fields.h"

namespace elements {
    template<typename T>
    using Number = BufferedComposite<fields::Number<T>>;

    template<typename T, size_t nind>
    struct NumberArray : BufferedComposite<fields::NumberArray<T, nind>> {
        NumberArray() : BufferedComposite<fields::NumberArray<T, nind>>(
                "Working number array") {}
    };

    struct Bitset : BufferedComposite<fields::Bitset> {
        Bitset(size_t nbit) : BufferedComposite<fields::Bitset>(
                nbit, "Working bitset") {}
    };

    struct FermionOnv : BufferedComposite<fields::FermionOnv> {
        FermionOnv(size_t nsite) : BufferedComposite<fields::FermionOnv>(
                nsite, "Working determinant") {}
    };

    struct BosonOnv : BufferedComposite<fields::BosonOnv> {
        BosonOnv(size_t nmode) : BufferedComposite<fields::BosonOnv>(
                nmode, "Working Boson ONV") {}
    };

    struct Onv : BufferedComposite<fields::Onv> {
        Onv(size_t nsite, size_t nmode) : BufferedComposite<fields::Onv>(
                nsite, nmode, "Working configuration") {}
    };
}

#endif //M7_ELEMENTS_H
#endif //M7_ELEMENTS_H
