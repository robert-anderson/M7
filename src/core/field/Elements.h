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
    struct FermionOnv : Element<fields::FermionOnv> {
        FermionOnv(size_t nsite):
        Element<fields::FermionOnv>(FermionOnvSpecifier(nsite), "Working determinant"){}
        using specs::FermionOnv::view_t::operator=;
    };

    struct BosonOnv : Element<fields::BosonOnv> {
        BosonOnv(size_t nmode):
                Element<fields::BosonOnv>(BosonOnvSpecifier(nmode), "Working boson ONV"){}
        using specs::BosonOnv::view_t::operator=;
    };

    struct FermiBosOnv : Element<fields::FermiBosOnv> {
        FermiBosOnv(size_t nsite, size_t nmode):
        Element<fields::FermiBosOnv>(nsite, nmode, "Working fermion-boson configuration"){}
    };

    using Onv = std::conditional<defs::bosons, FermiBosOnv, FermionOnv>::type;
}

#endif //M7_ELEMENTS_H
