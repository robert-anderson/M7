//
// Created by rja on 10/02/2021.
//

#ifndef M7_ELEMENTSZ_H
#define M7_ELEMENTSZ_H

#include "NdMultiFieldZ.h"
#include "BufferedTableZ.h"



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

template<typename ...Args>
struct ElementZ : WrappedRowZ, NdMultiFieldZ<0ul, Args...> {
    Buffer m_buffer;
    TableBaseZ m_table;
    ElementZ(Args &&... subfields) :
            NdMultiFieldZ<0ul, Args...>(&m_wrapped_row, {}, std::move(subfields)...),
            m_buffer("", 1), m_table(m_wrapped_row.m_dsize) {
        m_table.set_buffer(&m_buffer);
        m_wrapped_row.m_table_bw = &m_table.m_bw;
        m_wrapped_row.m_table_hwm = &m_table.m_hwm;
        m_table.push_back();
        m_wrapped_row.restart();
        NdMultiFieldZ<0ul, Args...>::restart();
    }
};



template<typename field_t>
struct NdItemZ : WrappedRowZ, NdFieldZ<0ul, field_t> {
    Buffer m_buffer;
    TableBaseZ m_table;
    NdItemZ(field_t&& field) : NdFieldZ<0ul, field_t>(&m_wrapped_row, {}, std::move(field)),
    m_buffer("", 1), m_table(m_wrapped_row.m_dsize) {
        m_table.set_buffer(&m_buffer);
        m_wrapped_row.m_table_bw = &m_table.m_bw;
        m_wrapped_row.m_table_hwm = &m_table.m_hwm;
        m_table.push_back();
        m_wrapped_row.restart();
    }
};


namespace itemsz {


    using FermionOnv = NdItemZ<FermionOnvFieldZ<0ul>>;

#if 0
    template<typename T, size_t nind>
    using number_array = ElementZ<NumberFieldZ<T, nind>>;

    struct FermionOnv : ElementZ<FermionOnvFieldZ> {
        FermionOnv(size_t nsite): ElementZ<FermionOnvFieldZ>(FermionOnvFieldZ(nsite)){}
    };

    struct FermiBosOnv : ElementZ<FermionOnvFieldZ, BosonOnvFieldZ> {
        FermiBosOnv(size_t nsite):
                ElementZ<FermionOnvFieldZ, BosonOnvFieldZ>(FermionOnvFieldZ(nsite), BosonOnvFieldZ(nsite)){}
    };
#endif
}


#endif //M7_ELEMENTSZ_H
