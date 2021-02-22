//
// Created by rja on 10/02/2021.
//

#ifndef M7_BUFFEREDFIELDS_H
#define M7_BUFFEREDFIELDS_H

#include "src/core/field/MultiField.h"
#include "BufferedTable.h"
#include "src/core/field/Fields.h"


//template<typename row_t>
//struct RowWrapperZ {
//    row_t m_row;
//
//};

/*
template<size_t nind, typename ...Args>
struct SingleFieldRowZ : Row {
    NdMultiFieldZ<nind, ...Args> m_field;
    SingleFieldRowZ(std::array<size_t, nind> shape, Args&&... subfields) :
    Row(), m_field(this, shape, subfields...) {}
};

template<typename nd_field_t>
struct BufferedSingleFieldTableZ : BufferedTable<SingleFieldRowZ<nd_field_t>> {
    template<typename ...Args>
    SingleFieldTableZ(Args... args):
            Table<SingleFieldRowZ<nd_field_t>>(SingleFieldRowZ<nd_field_t>(args...)){}
};

template<typename nd_field_t>
struct ElementZ : nd_field_t, BufferedTable<row_t> {
    ElementZ(Args... args): NdFieldZ<0ul, row_t>()
};
*/

struct WrappedRowZ {
    Row m_wrapped_row;
};

template<typename multi_field_t>
struct BufferedMultiFieldz : WrappedRowZ, multi_field_t {
    Buffer m_buffer;
    TableBase m_table;
    BufferedMultiFieldz(multi_field_t multi_field) :
            multi_field_t(multi_field),
            m_buffer("", 1),
            m_table((multi_field_t::add_to_row(&m_wrapped_row), m_wrapped_row.m_dsize)) {
        m_table.set_buffer(&m_buffer);
        m_wrapped_row.m_table_bw = &m_table.m_bw;
        m_wrapped_row.m_table_hwm = &m_table.m_hwm;
        m_table.push_back();
        m_wrapped_row.restart();
    }
};

template<typename field_t>
struct BufferedFieldZ : WrappedRowZ, field_t {
    static_assert(std::is_base_of<FieldBase, field_t>::value, "Template arg must be derived from FieldBase");
    Buffer m_buffer;
    TableBase m_table;
    BufferedFieldZ(field_t field) : field_t(field),
        m_buffer("", 1),
        m_table((field_t::add_to_row(&m_wrapped_row), m_wrapped_row.m_dsize)) {
        m_table.set_buffer(&m_buffer);
        m_wrapped_row.m_table_bw = &m_table.m_bw;
        m_wrapped_row.m_table_hwm = &m_table.m_hwm;
        m_table.push_back();
        m_wrapped_row.restart();
    }
};


namespace buffered {


    struct FermionOnv : BufferedFieldZ<fields::FermionOnv> {
        using fields::FermionOnv::operator=;
        FermionOnv(size_t nsite) : BufferedFieldZ<fields::FermionOnv>({nullptr, nsite}){}
    };

    struct FermionOnvs : BufferedFieldZ<fields::FermionOnvs> {
        FermionOnvs(size_t nitem, size_t nsite) : BufferedFieldZ<fields::FermionOnvs>({nullptr, nitem, nsite}){}
    };

    struct BosonOnv : BufferedFieldZ<fields::BosonOnv> {
        using fields::BosonOnv::operator=;
        BosonOnv(size_t nsite) : BufferedFieldZ<fields::BosonOnv>({nullptr, nsite}){}
    };

    struct BosonOnvs : BufferedFieldZ<fields::BosonOnvs> {
        BosonOnvs(size_t nitem, size_t nsite) : BufferedFieldZ<fields::BosonOnvs>({nullptr, nitem, nsite}){}
    };

    struct FermiBosOnv : BufferedMultiFieldz<fields::FermiBosOnv> {
        using fields::FermiBosOnv::operator=;
        FermiBosOnv(size_t nsite):
        BufferedMultiFieldz<fields::FermiBosOnv>({nullptr, nsite}){}
    };
    using FermiBosOnvs = BufferedMultiFieldz<fields::FermiBosOnvs>;

    template<bool enable_bosons=defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;

    template<bool enable_bosons=defs::enable_bosons>
    using Onvs = typename std::conditional<enable_bosons, FermiBosOnvs, FermionOnvs>::type;

}


#endif //M7_BUFFEREDFIELDS_H
