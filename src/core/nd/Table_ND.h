//
// Created by rja on 21/10/2020.
//

#ifndef M7_TABLE_ND_H
#define M7_TABLE_ND_H


#include <iostream>
#include <array>
#include <vector>
#include <climits>
#include <cstring>
#include <src/core/util/utils.h>
#include <src/core/util/defs.h>
#include "NdFormat.h"

struct Field;

struct Table {
    size_t m_row_size;
    std::vector<Field *> m_fields;
    char* m_data;

    char *begin() {
        return m_data;
    }

    char *begin(const size_t &irow) {
        return begin() + irow * m_row_size;
    }

    size_t add_field(const Field& field){
        return 0;// offset
    }
};

struct Field {
    Table *m_table;
    const size_t m_nelement;
    const size_t m_element_size;
    const size_t m_size;
    const size_t m_offset;
    const std::string m_description;

    Field(Table *table, size_t nelement, size_t element_size, std::string description) :
            m_table(table), m_nelement(nelement), m_element_size(element_size), m_size(nelement * element_size),
            m_offset(table->add_field(*this)), m_description(std::move(description)) {}

    char *begin(const size_t &irow) const {
        return m_table->begin(irow);
    }

    char *raw_ptr(const size_t &irow, const size_t &iflat) const {
        return begin(irow) + iflat * m_element_size;
    }

    template<typename ...Args>
    std::pair<const char *, size_t> raw_view(const size_t &irow, Args...inds) const {
        return {raw_ptr(irow, inds...), m_element_size};
    }
};


template<typename field_t>
struct RowAccessor {
    const field_t *m_field;
    char *m_ptr;

    RowAccessor(field_t *field, size_t irow) : m_field(field), m_ptr(field->begin(irow)) {}

    template<typename ...Args>
    typename field_t::view_t operator()(Args... inds) {
        return field_t::view_factory(m_ptr);
    }
};


template<size_t nind>
struct NdField : Field {
    NdFormat<nind> m_format;

    NdField(Table *table, std::array<size_t, nind> shape, size_t element_size, std::string description) :
            Field(table, NdFormat<nind>(shape).nelement(), element_size, description),
            m_format(shape) {}

    template<typename ...Args>
    char *raw_ptr(const size_t &irow, Args...inds) const {
        return Field::raw_ptr(irow, m_format.flat(inds...));
    }

    template<typename ...Args>
    std::pair<char *, size_t> raw_view(const size_t &irow, Args...inds) {
        return {raw_ptr(irow, inds...), m_element_size};
    }
};

template<size_t nind, size_t nind_view>
struct NdArrayField : NdField<nind> {
    static_assert(nind_view, "If the view is scalar, use NdField instead");
    NdFormat<nind_view> m_view_format;

    NdArrayField(Table *table, std::array<size_t, nind> shape,
                 std::array<size_t, nind_view> view_shape, size_t element_size, std::string description) :
            NdField<nind>(table, shape, NdFormat<nind_view>(view_shape).nelement() * element_size, description),
            m_view_format(view_shape) {}

};


template<typename T, size_t nind>
struct NumericField : NdField<nind> {
    NumericField(Table *table, std::array<size_t, nind> shape, std::string description) :
            NdField<nind>(table, shape, sizeof(T), description) {}

    template<typename ...Args>
    T &operator()(const size_t &irow, Args... inds) {
        return (T *) raw_ptr(irow, inds...);
    }
};

template<typename T, size_t nind, size_t nind_view>
struct NumericArraySpecifier : NdArrayField<nind, nind_view> {
    NumericArraySpecifier(Table *table, std::array<size_t, nind> shape, std::array<size_t, nind_view> view_shape, std::string description) :
            NdArrayField<nind, nind_view>(table, shape, view_shape, sizeof(T), description) {}

    struct NdView {
        const NumericArraySpecifier& m_field;
        char* m_ptr;
        NdView(const NumericArraySpecifier& field, const size_t &irow, const size_t &iflat):
                m_field(field), m_ptr(field.raw_ptr(irow, iflat)){}

        template<typename ...Args>
        T& operator()(Args... inds){
            return ((T*)m_ptr)[m_field.m_view_format.flat(inds...)];
        }

        NdView& operator=(const NdView& other){
            if (&other!=this)
                std::memcpy(m_ptr, other.m_ptr, m_field.m_element_size);
            return *this;
        }
    };

    template<typename ...Args>
    NdView operator()(const size_t &irow, Args... inds) {
        return NdView(*this, irow, NdField<nind>::m_format.flat(inds...));
    }
};

template<size_t nind>
struct BitsetField : NdField<nind> {
    const size_t m_nbit;
    BitsetField(Table *table, std::array<size_t, nind> shape, size_t nbit, std::string description) :
            NdField<nind>(table, shape, integer_utils::divceil(nbit, defs::nbit_data), description),
                    m_nbit(nbit){}
    struct View {
        const BitsetField& m_field;
        defs::data_t* m_ptr;
        View(const BitsetField& field, char* ptr): m_field(field), m_ptr((defs::data_t*)ptr){}

    };
//
//    template<typename ...Args>
//    View operator()(const size_t &irow, Args... inds) {
//        return View(*this, irow, NdField<nind>::m_format.flat(inds...));
//    }
//    struct NdView {
//        const BitsetField& m_field;
//        char* m_ptr;
//        NdView(const BitsetField& field, const size_t &irow, const size_t &iflat):
//                m_field(field), m_ptr(field.raw_ptr(irow, iflat)){}
//
//        template<typename ...Args>
//        View operator()(Args... inds){
//            return View(m_ptr+m_field.m_view_format.flat(inds...)*m_field.m_element_size);
//        }
//
//        NdView& operator=(const NdView& other){
//            if (&other!=this)
//                std::memcpy(m_ptr, other.m_ptr, m_field.m_element_size);
//            return *this;
//        }
//    };
//
//    template<typename ...Args>
//    NdView operator()(const size_t &irow, Args... inds) {
//        return NdView(*this, irow, NdField<nind>::m_format.flat(inds...));
//    }
};






#endif //M7_TABLE_ND_H
