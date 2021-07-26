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
            m_table((multi_field_t::add_to_row(&m_wrapped_row), m_wrapped_row.m_size)) {
        m_table.set_buffer(&m_buffer);
        m_wrapped_row.m_table = &m_table;
        m_table.push_back();
        m_wrapped_row.restart();
    }
};

template<typename field_t>
struct BufferedField : WrappedRow, field_t {
    static_assert(std::is_base_of<FieldBase, field_t>::value, "Template arg must be derived from FieldBase");
    Buffer m_buffer;
    TableBase m_table;
    BufferedField(const field_t& field) : field_t(field),
                                   m_buffer("", 1),
                                   m_table((field_t::add_to_row(&m_wrapped_row), m_wrapped_row.m_size)) {
        m_table.set_buffer(&m_buffer);
        m_wrapped_row.m_table = &m_table;
        m_table.push_back();
        m_wrapped_row.restart();
    }

    BufferedField(const BufferedField& other): BufferedField(static_cast<const field_t&>(other)){
        // copy data
        *this = other;
    }

    BufferedField& operator=(const field_t& other){
        field_t::operator=(other);
        return *this;
    }

    BufferedField& operator=(const BufferedField& other){
        field_t::operator=(other);
        return *this;
    }
};


namespace buffered {

    template<typename T, size_t nind>
    struct NdBitset : BufferedField<fields::NdBitset<T, nind>> {
        typedef BufferedField<fields::NdBitset<T, nind>> base_t;
        typedef typename fields::NdBitset<T, nind>::inds_t inds_t;
        using fields::NdBitset<T, nind>::operator=;
        NdBitset(inds_t shape) : BufferedField<fields::NdBitset<T, nind>>({nullptr, shape}){}
        NdBitset(const fields::NdBitset<T, nind>& field): BufferedField<fields::NdBitset<T, nind>>(field){}
    };


    template<typename T>
    struct Bitset : NdBitset<T, 1ul> {
        Bitset(size_t nbit): NdBitset<T, 1ul>({nbit}){}
    };


    template<typename T, size_t nind>
    struct Numbers : BufferedField<fields::Numbers<T, nind>> {
        typedef BufferedField<fields::Numbers<T, nind>> base_t;
        typedef typename fields::Numbers<T, nind>::inds_t inds_t;
        using fields::Numbers<T, nind>::operator=;
        Numbers(inds_t shape) : base_t({nullptr, shape}){}
        Numbers(inds_t shape, T init_value) : base_t({nullptr, shape}){
            *this = init_value;
        }
        Numbers(const fields::Numbers<T, nind>& field): BufferedField<fields::Numbers<T, nind>>(field){}
        Numbers(const Numbers& field): BufferedField<fields::Numbers<T, nind>>(field){}
        Numbers& operator=(const Numbers& field){
            base_t::operator=(field);
            return *this;
        }
    };

    template<typename T>
    struct Number : BufferedField<fields::Number<T>> {
        using BufferedField<fields::Number<T>>::operator=;
        Number(): BufferedField<fields::Number<T>>({{}}){}
        operator T&(){return (*this)[0];}
        operator const T&() const {return (*this)[0];}
        Number& operator=(const T& v){
            static_cast<T&>(*this) = v;
            return *this;
        }
    };

    struct FrmOnv : BufferedField<fields::FrmOnv> {
        using fields::FrmOnv::operator=;
        FrmOnv(size_t nsite) : BufferedField<fields::FrmOnv>({nullptr, nsite}){}
        FrmOnv& operator=(const FrmOnv& other){
            base_t::operator=(other);
            return *this;
        }
        FrmOnv(const FrmOnv& other): FrmOnv(other.m_nsite){}
    };

    struct BosOnv : BufferedField<fields::BosOnv> {
        using fields::BosOnv::operator=;
        BosOnv(size_t nsite) : BufferedField<fields::BosOnv>({nullptr, nsite}){}
    };


    struct FrmBosOnv : BufferedMultiField<fields::FrmBosOnv> {
        using fields::FrmBosOnv::operator=;
        FrmBosOnv(size_t nsite):
                BufferedMultiField<fields::FrmBosOnv>({nullptr, nsite}){}
    };

    template<bool enable_bosons=defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FrmBosOnv, FrmOnv>::type;

    typedef std::tuple<FrmOnv, FrmBosOnv> mbf_tup_t;
    typedef typename std::tuple_element<defs::mbf_type_id, mbf_tup_t>::type mbf_t;

    struct FermionMevInds : BufferedMultiField<fields::FermionMevInds> {
        using fields::FermionMevInds::operator=;
        using fields::FermionMevInds::m_ann;
        using fields::FermionMevInds::m_cre;
        FermionMevInds(size_t nann, size_t ncre):
                BufferedMultiField<fields::FermionMevInds>({nullptr, nann, ncre}){}

        FermionMevInds(size_t nop): FermionMevInds(nop, nop){}
    };
}

template<typename field_t>
struct SingleFieldRow : Row {
    field_t m_field;
    template<typename ...Args>
    SingleFieldRow(Args... args): Row(), m_field(this, args...){}

    field_t &key_field() {
        return m_field;
    };

};


#endif //M7_BUFFEREDFIELDS_H
