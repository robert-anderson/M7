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
struct SingleFieldTable: Table {
    static_assert(std::is_base_of<NdFieldGroup<0ul>, viewable_t>::value, "Template arg must be a scalar NdFieldGroup");
    viewable_t m_field;
	template <typename ...Args>
    SingleFieldTable(Args... args):
    Table(),
    m_field(this, args..., "Working field"){}
};

template<typename viewable_t>
struct BufferedSingleFieldTable: BufferedTable<SingleFieldTable<viewable_t>>{
    typedef BufferedTable<SingleFieldTable<viewable_t>> base_t;
	template <typename ...Args>
    BufferedSingleFieldTable(Args... args): base_t("", args...){
        base_t::resize(1ul);
        base_t::push_back();
    }
};

template <typename viewable_t>
struct Element : BufferedSingleFieldTable<viewable_t>, viewable_t::view_t {
    typedef BufferedSingleFieldTable<viewable_t> base_t;
	template <typename ...Args>
    Element(Args... args):
    base_t(args...), viewable_t::view_t(base_t::m_field(0)){}
};

namespace elements {
    struct FermionOnv : Element<fields::Det> {
        FermionOnv(size_t nsite):
            Element<fields::Det>({nsite}){
        }
        FermionOnv(const views::Det & view): FermionOnv(view.nsite()){
            *this = view;
        }
        using specs::FermionOnv::view_t::operator=;
    };

    struct BosonOnv : Element<fields::BosonOnv> {
        BosonOnv(size_t nmode):
            Element<fields::BosonOnv>({nmode}){}
        BosonOnv(const views::BosonOnv& view): BosonOnv(view.nmode()){
            *this = view;
        }
        using specs::BosonOnv::view_t::operator=;
    };

    struct FermiBosOnv : Element<fields::Onv<1>> {
        FermiBosOnv(size_t nsite):
            Element<fields::Onv<1>>(nsite){}
        FermiBosOnv(const views::Onv<1>& view):
        FermiBosOnv(view.m_fonv.nsite()){
            *this = view;
        }
        using fields::Onv<1>::view_t::operator=;
    };

    template<bool enable_bosons=defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;

    using Det = Onv<0>;
}

#endif //M7_ELEMENTS_H
