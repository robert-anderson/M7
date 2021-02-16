//
// Created by rja on 09/02/2021.
//

#ifndef M7_NDMULTIFIELDZ_H
#define M7_NDMULTIFIELDZ_H

#include "RowZ.h"
#include "src/core/nd/NdSelector.h"

struct RowFieldBaseZ {
    RowZ *m_row;
    RowFieldBaseZ(RowZ* row): m_row(row){}

    RowFieldBaseZ(const RowFieldBaseZ& other): m_row(other.m_row ? other.m_row->m_child : nullptr){}
};

template<size_t nind_item, typename ...Args>
struct NdMultiFieldZ : RowFieldBaseZ {
    NdFormat<nind_item> m_format;
    std::tuple<Args...> m_subfields;
    std::vector<char> m_null_field_string;

    NdMultiFieldZ(RowZ *row, std::array<size_t, nind_item> shape, Args&&... subfields) :
            RowFieldBaseZ(row), m_format(shape), m_subfields(subfields...) {
        init();
        m_null_field_string.assign(max_size(), 0);
    }

    /*
     * NdMultiFields must be bound to rows, and so all copies are assumed to be done as
     * part of a row copy. This is in contrast to the FieldBase class, which may be
     * legitimately copied without reference to a row.
     */
    NdMultiFieldZ(const NdMultiFieldZ& other) :
            RowFieldBaseZ(other.m_row->m_child),
            m_format(other.m_format),
            m_subfields(other.m_subfields){
        init();
    }

    NdMultiFieldZ& operator=(const NdMultiFieldZ& other) {
        MPI_ASSERT(is_comparable(other), "Shapes are incompatible");
        struct fn_t {
            const NdMultiFieldZ& other;
            void operator()(FieldBaseZ &f1, const FieldBaseZ &f2) { f1=f2; }
        };
        fn_t fn;
        tuple_utils::for_each_pair(m_subfields, other.m_subfields, fn);
        return *this;
    }

    void init() {
        if (!m_row) return;
        struct fn_t {
            RowZ* m_row;
            NdFormat<nind_item>* m_format;
            void operator()(ItemFormattedFieldBaseZ<nind_item> &f) {
                f.m_item_format = m_format;
                f.m_nitem = m_format->nelement();
                f.m_size = f.m_item_size*f.m_nitem;
                f.m_row_offset = m_row->add_field(&f);
                f.m_row = m_row;
            }
        };
        fn_t fn {m_row, &m_format};
        tuple_utils::for_each(m_subfields, fn);
    }

    bool operator==(const NdMultiFieldZ &other) const {
        struct fn_t {
            bool m_and = true;

            void operator()(const FieldBaseZ &f1, const FieldBaseZ &f2) { m_and &= f1==f2; }
        };
        fn_t fn;
        tuple_utils::for_each_pair(m_subfields, other.m_subfields, fn);
        return fn.m_and;
    }

    bool is_comparable(const NdMultiFieldZ &other) const {
        struct fn_t {
            bool m_and = true;
            void operator()(const FieldBaseZ &f1, const FieldBaseZ &f2) { m_and &= f1.is_comparable(f2); }
        };
        fn_t fn;
        tuple_utils::for_each_pair(m_subfields, other.m_subfields, fn);
        return fn.m_and;
    }



protected:
    bool max_size() const {
        struct fn_t {
            size_t m_max = 0ul;
            void operator()(const FieldBaseZ &f) { if (f.m_size>m_max) m_max = f.m_size; }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
        return fn.m_max;
    }

    bool is_zero() const {
        struct fn_t {
            bool m_and = true;
            void operator()(const FieldBaseZ &f) { m_and &= f.is_zero(); }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
        return fn.m_and;
    }

    defs::hash_t hash() const {
        struct fn_t {
            defs::hash_t m_hash = 0;
            void operator()(const FieldBaseZ &f) { m_hash ^= f.hash(); }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
        return fn.m_hash;
    }

    std::string to_string() const {
        struct fn_t {
            std::string m_str;

            void operator()(const FieldBaseZ &f) { m_str += f.to_string() + " "; }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
        return fn.m_str;
    }

    template<size_t ifield>
    const typename std::tuple_element<ifield, std::tuple<Args...>>::type &get() const {
        return std::get<ifield>(m_subfields);
    }

    template<size_t ifield>
    typename std::tuple_element<ifield, std::tuple<Args...>>::type &get() {
        return std::get<ifield>(m_subfields);
    }
};


template<size_t nind_item, typename field_t>
struct NdFieldZ : RowFieldBaseZ, field_t {
    using RowFieldBaseZ::m_row;
    using field_t::operator=;
    NdFormat<nind_item> m_format;
    NdFieldZ(RowZ *row, std::array<size_t, nind_item> item_shape, field_t&& field) :
            RowFieldBaseZ(row), field_t(field), m_format(item_shape){
        static_cast<ItemFormattedFieldBaseZ<nind_item>&>(*this).m_item_format = &m_format;
        FieldBaseZ::m_nitem = m_format.nelement();
        FieldBaseZ::m_size = this->m_item_size*this->m_nitem;
        FieldBaseZ::m_row_offset = m_row->add_field(this);
        FieldBaseZ::m_row = m_row;
    }

    NdFieldZ(const NdFieldZ& other):
            RowFieldBaseZ(other.m_row->m_child), field_t(other), m_format(other.m_format.shape()){
        static_cast<ItemFormattedFieldBaseZ<nind_item>&>(*this).m_item_format = &m_format;
        FieldBaseZ::m_nitem = m_format.nelement();
        FieldBaseZ::m_size = this->m_item_size*this->m_nitem;
        ASSERT(m_row);
        FieldBaseZ::m_row_offset = m_row->add_field(this);
        FieldBaseZ::m_row = m_row;
    }

};


#endif //M7_NDMULTIFIELDZ_H