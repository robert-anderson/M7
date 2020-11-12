//
// Created by rja on 01/11/2020.
//

#ifndef M7_ELEMENTS_H
#define M7_ELEMENTS_H

#include "src/core/table/Table.h"
#include "src/core/table/BufferedTable.h"
#include "Fields.h"
#include "Views.h"

template<typename viewable_t>
struct SingleFieldTable: TableX {
    static_assert(std::is_base_of<NdFieldGroup<0ul>, viewable_t>::value, "Template arg must be a scalar NdFieldGroup");
    typedef typename viewable_t::params_t params_t;
    viewable_t m_field;
    SingleFieldTable(params_t p, std::string description):
    TableX(),
    m_field(this, p, description){}
};

template<typename viewable_t>
struct BufferedSingleFieldTable: BufferedTable<SingleFieldTable<viewable_t>>{
    typedef BufferedTable<SingleFieldTable<viewable_t>> base_t;
    typedef typename viewable_t::params_t params_t;
    BufferedSingleFieldTable(params_t p, std::string description): base_t(p, description){
        base_t::expand(1ul);
        base_t::push_back();
    }
};

template <typename viewable_t>
struct Element : BufferedSingleFieldTable<viewable_t>, viewable_t::view_t {
    typedef BufferedSingleFieldTable<viewable_t> base_t;
    typedef typename viewable_t::params_t params_t;
    Element(params_t p, std::string description):
    base_t(p, description),
    viewable_t::view_t(base_t::m_field(0)){}
};

namespace elements {
    struct FermionOnv : Element<fields::FermionOnv> {
        FermionOnv(size_t nsite):
            Element<fields::FermionOnv>({nsite}, "Working determinant"){}
        FermionOnv(fields::FermionOnv::params_t p): FermionOnv(p.m_nsite){}
        FermionOnv(const views::FermionOnv& view): FermionOnv(view.nsite()){
            *this = view;
        }
        using specs::FermionOnv::view_t::operator=;
    };

    struct BosonOnv : Element<fields::BosonOnv> {
        BosonOnv(size_t nmode):
            Element<fields::BosonOnv>({nmode}, "Working boson ONV"){}
        BosonOnv(fields::BosonOnv::params_t p): BosonOnv(p.m_nmode){}
        BosonOnv(const views::BosonOnv& view): BosonOnv(view.nmode()){
            *this = view;
        }
        using specs::BosonOnv::view_t::operator=;
    };

    struct FermiBosOnv : Element<fields::FermiBosOnv> {
        FermiBosOnv(size_t nsite, size_t nmode):
            Element<fields::FermiBosOnv>({nsite, nmode}, "Working fermion-boson configuration"){}
        FermiBosOnv(fields::FermiBosOnv::params_t p): FermiBosOnv(p.m_nsite, p.m_nmode){}
        FermiBosOnv(const views::FermiBosOnv& view):
        FermiBosOnv(view.m_fonv.nsite(), view.m_bonv.nmode()){
            *this = view;
        }
        using fields::FermiBosOnv::view_t::operator=;
    };

    using Onv = std::conditional<defs::bosons, FermiBosOnv, FermionOnv>::type;
}

#endif //M7_ELEMENTS_H
