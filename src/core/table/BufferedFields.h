//
// Created by rja on 10/02/2021.
//

#ifndef M7_BUFFEREDFIELDS_H
#define M7_BUFFEREDFIELDS_H

#include "src/core/field/MultiField.h"
#include "BufferedTable.h"
#include "src/core/field/Fields.h"

struct WrappedRow {
    Row m_wrapped_row;
};

template<typename multi_field_t>
struct BufferedMultiField : WrappedRow, multi_field_t {
    Buffer m_buffer;
    TableBase m_table;
    BufferedMultiField(multi_field_t multi_field) :
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
struct BufferedField : WrappedRow, field_t {
    static_assert(std::is_base_of<FieldBase, field_t>::value, "Template arg must be derived from FieldBase");
    Buffer m_buffer;
    TableBase m_table;
    BufferedField(field_t field) : field_t(field),
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

    template<typename T>
    struct Vector : BufferedField<fields::Vector<T>> {
        Vector(size_t nelement) : BufferedField<fields::Vector<T>>({nullptr, nelement}){}
    };

    struct FermionOnv : BufferedField<fields::FermionOnv> {
        using fields::FermionOnv::operator=;
        FermionOnv(size_t nsite) : BufferedField<fields::FermionOnv>({nullptr, nsite}){}
        FermionOnv& operator=(const FermionOnv& other){
            fields::FermionOnv::operator=(other);
            return *this;
        }
    };

    struct FermionOnvs : BufferedField<fields::FermionOnvs> {
        FermionOnvs(size_t nitem, size_t nsite) : BufferedField<fields::FermionOnvs>({nullptr, nitem, nsite}){}
    };

    struct BosonOnv : BufferedField<fields::BosonOnv> {
        using fields::BosonOnv::operator=;
        BosonOnv(size_t nsite) : BufferedField<fields::BosonOnv>({nullptr, nsite}){}
    };

    struct BosonOnvs : BufferedField<fields::BosonOnvs> {
        BosonOnvs(size_t nitem, size_t nsite) : BufferedField<fields::BosonOnvs>({nullptr, nitem, nsite}){}
    };

    struct FermiBosOnv : BufferedMultiField<fields::FermiBosOnv> {
        using fields::FermiBosOnv::operator=;
        FermiBosOnv(size_t nsite):
                BufferedMultiField<fields::FermiBosOnv>({nullptr, nsite}){}
    };
    using FermiBosOnvs = BufferedMultiField<fields::FermiBosOnvs>;

    template<bool enable_bosons=defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;

    template<bool enable_bosons=defs::enable_bosons>
    using Onvs = typename std::conditional<enable_bosons, FermiBosOnvs, FermionOnvs>::type;

}


#endif //M7_BUFFEREDFIELDS_H
