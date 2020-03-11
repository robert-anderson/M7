//
// Created by Robert John Anderson on 2020-03-09.
//

#ifndef M7_FIELD_H
#define M7_FIELD_H

#include <src/multidim/Indexer.h>
#include "FieldSet.h"
#include <cstring>


template<typename T, size_t nind = 1>
struct Field : public FieldSet::FieldBase {
    const Indexer<nind> m_indexer;

    template<typename ...Args>
    Field(FieldSet *field_set, Args... shape) :
            FieldBase(field_set, FieldSet::nbit_per_element<T>(), Indexer<nind>(shape...).nelement(), typeid(T)),
            m_indexer(Indexer<nind>(shape...)) {}

private:
    inline T *flat_get(const size_t &irow, const size_t &flat) {
        assert(flat < m_nelement);
        return (T *) (m_field_set->m_buffer + irow * m_field_set->m_length +
                      m_offset.first) + m_offset.second + flat;
    }

public:
    template<typename U=T, typename ...Args>
    typename std::enable_if<!std::is_same<U, void>::value, U *>::type
    zero(const size_t &irow){
        std::memset(flat_get(irow, 0), 0, m_nelement*sizeof(U));
    }


    template<typename U=T, typename ...Args>
    typename std::enable_if<!std::is_same<U, bool>::value, U *>::type
    operator()(const size_t &irow, Args... inds) {
        return flat_get(irow, m_indexer.get(inds...));
    }

    template<typename ...Args>
    std::string to_string(size_t irow, Args... inds) {
        return std::to_string(*(*this)(irow, inds...));
    }

    std::string to_string(size_t irow) override {
        std::string out = "";
        for (size_t i = 0ul; i < m_nelement; ++i) out += std::to_string(*flat_get(irow, i))+" ";
        return out;
    }
};

template<typename T>
static inline void clr_bit(T &x, size_t i) {
    x &= ~((T) 1ul << i);
}

template<typename T>
static inline void set_bit(T &x, size_t i) {
    x |= ((T) 1ul << i);
}

template<typename T>
static inline bool get_bit(T &x, size_t i) {
    return (x >> i) & T(1ul);
}

template<size_t nind = 1>
struct Flag : public Field<bool, nind> {
    template<typename ...Args>
    Flag(FieldSet *field_set, Args... shape) : Field<bool, nind>(field_set, shape...) {}

    using Field<bool, nind>::m_field_set;
    using Field<bool, nind>::m_offset;

    template<typename ...Args>
    void set(const size_t &irow, Args... inds) {
        set_bit(*(m_field_set->m_buffer + irow * m_field_set->m_length + m_offset.first), m_offset.second);
    }

    template<typename ...Args>
    void clr(const size_t &irow, Args... inds) {
        clr_bit(*(m_field_set->m_buffer + irow * m_field_set->m_length + m_offset.first), m_offset.second);
    }

    template<typename ...Args>
    bool get(const size_t &irow, Args... inds) {
        return get_bit(*(m_field_set->m_buffer + irow * m_field_set->m_length + m_offset.first), m_offset.second);
    }
};


template<size_t nind = 1>
struct Bitfield : public Field<defs::data_t, nind+1> {
    const size_t m_nbit;

    template<typename ...Args>
    Bitfield(FieldSet *field_set, size_t nbit, Args... shape) : Field<defs::data_t, nind+1>(
            field_set, shape..., integer_utils::divceil(nbit, FieldSet::nbit_per_element<defs::data_t>())),
            m_nbit(nbit) {}

    using Field<defs::data_t, nind+1>::m_field_set;
    using Field<defs::data_t, nind+1>::m_offset;

    template<typename ...Args>
    void set(const size_t &irow, const size_t &ibit, Args... inds) {
        set_bit(*(*this)(inds..., ibit/FieldSet::nbit_per_element<defs::data_t>()),
            ibit%FieldSet::nbit_per_element<defs::data_t>());
    }

    template<typename ...Args>
    void clr(const size_t &irow, const size_t &ibit, Args... inds) {
        clr_bit(*(*this)(inds..., ibit/FieldSet::nbit_per_element<defs::data_t>()),
                ibit%FieldSet::nbit_per_element<defs::data_t>());
    }

    template<typename ...Args>
    bool get(const size_t &irow, const size_t &ibit, Args... inds) {
        return get_bit(*(*this)(inds..., ibit/FieldSet::nbit_per_element<defs::data_t>()),
                ibit%FieldSet::nbit_per_element<defs::data_t>());
    }

    /*
    template<typename ...Args>
    std::string to_string(size_t irow, Args... inds) override {
        return std::to_string(*(*this)(irow, inds...));
    }*/
};


template<size_t nind = 1>
struct DeterminantField : public Bitfield<nind+1>{
    template<typename ...Args>
    DeterminantField(FieldSet *field_set, size_t nbit, Args... shape) : Bitfield<nind+1>(
            field_set, nbit, shape..., 2){}
};


#endif //M7_FIELD_H
