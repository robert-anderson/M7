//
// Created by rja on 10/02/2021.
//

#ifndef M7_BUFFEREDFIELDS_H
#define M7_BUFFEREDFIELDS_H

#include "MultiFieldZ.h"
#include "BufferedTableZ.h"
#include "FieldsZ.h"


//template<typename row_t>
//struct RowWrapperZ {
//    row_t m_row;
//
//};

/*
template<size_t nind, typename ...Args>
struct SingleFieldRowZ : RowZ {
    NdMultiFieldZ<nind, ...Args> m_field;
    SingleFieldRowZ(std::array<size_t, nind> shape, Args&&... subfields) :
    RowZ(), m_field(this, shape, subfields...) {}
};

template<typename nd_field_t>
struct BufferedSingleFieldTableZ : BufferedTableZ<SingleFieldRowZ<nd_field_t>> {
    template<typename ...Args>
    SingleFieldTableZ(Args... args):
            TableZ<SingleFieldRowZ<nd_field_t>>(SingleFieldRowZ<nd_field_t>(args...)){}
};

template<typename nd_field_t>
struct ElementZ : nd_field_t, BufferedTableZ<row_t> {
    ElementZ(Args... args): NdFieldZ<0ul, row_t>()
};
*/

struct WrappedRowZ {
    RowZ m_wrapped_row;
};

template<typename multi_field_t>
struct BufferedMultiFieldz : WrappedRowZ, multi_field_t {
    Buffer m_buffer;
    TableBaseZ m_table;
    BufferedMultiFieldz(multi_field_t&& multi_field) :
            multi_field_t(std::move(multi_field)),
            m_buffer("", 1), m_table(m_wrapped_row.m_dsize) {
        MPI_ASSERT(!multi_field_t::m_row, "MultiField must not be already associated with a row");
        multi_field_t::add_to_row(&m_wrapped_row);
        m_table.set_buffer(&m_buffer);
        m_wrapped_row.m_table_bw = &m_table.m_bw;
        m_wrapped_row.m_table_hwm = &m_table.m_hwm;
        m_table.push_back();
        m_wrapped_row.restart();
    }
};

template<typename field_t>
struct BufferedFieldZ : WrappedRowZ, field_t {
    static_assert(std::is_base_of<FieldBaseZ, field_t>::value, "Template arg must be derived from FieldBase");
    Buffer m_buffer;
    TableBaseZ m_table;
    BufferedFieldZ(field_t&& field) : field_t(std::move(field)),
        m_buffer("", 1), m_table(m_wrapped_row.m_dsize) {
        MPI_ASSERT(!FieldBaseZ::m_row, "Field must not be already associated with a row");
        FieldBaseZ::add_to_row(&m_wrapped_row);
        m_table.set_buffer(&m_buffer);
        m_wrapped_row.m_table_bw = &m_table.m_bw;
        m_wrapped_row.m_table_hwm = &m_table.m_hwm;
        m_table.push_back();
        m_wrapped_row.restart();
    }
};


namespace buffered {


    struct FermionOnv : BufferedFieldZ<fieldsz::FermionOnv> {
        using fieldsz::FermionOnv::operator=;
        FermionOnv(size_t nsite) : BufferedFieldZ<fieldsz::FermionOnv>({nullptr, nsite}){}
    };

    struct FermionOnvs : BufferedFieldZ<fieldsz::FermionOnvs> {
        FermionOnvs(size_t nitem, size_t nsite) : BufferedFieldZ<fieldsz::FermionOnvs>({nullptr, nitem, nsite}){}
    };

    struct BosonOnv : BufferedFieldZ<fieldsz::BosonOnv> {
        using fieldsz::BosonOnv::operator=;
        BosonOnv(size_t nsite) : BufferedFieldZ<fieldsz::BosonOnv>({nullptr, nsite}){}
    };

    struct BosonOnvs : BufferedFieldZ<fieldsz::BosonOnvs> {
        BosonOnvs(size_t nitem, size_t nsite) : BufferedFieldZ<fieldsz::BosonOnvs>({nullptr, nitem, nsite}){}
    };

    struct FermiBosOnv : BufferedMultiFieldz<fieldsz::FermiBosOnv> {
        using fieldsz::FermiBosOnv::operator=;
        FermiBosOnv(size_t nsite):
        BufferedMultiFieldz<fieldsz::FermiBosOnv>({nullptr, nsite}){}
    };
    using FermiBosOnvs = BufferedMultiFieldz<fieldsz::FermiBosOnvs>;

    template<bool enable_bosons=defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;

    template<bool enable_bosons=defs::enable_bosons>
    using Onvs = typename std::conditional<enable_bosons, FermiBosOnvs, FermionOnvs>::type;

}


#endif //M7_BUFFEREDFIELDS_H
