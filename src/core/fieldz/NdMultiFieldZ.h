//
// Created by rja on 09/02/2021.
//

#ifndef M7_NDMULTIFIELDZ_H
#define M7_NDMULTIFIELDZ_H

#include "RowZ.h"
#include "FieldBaseZ.h"
#include "src/core/nd/NdSelector.h"


template<typename ...Args>
struct MultiFieldZ {
    RowZ *m_row;
    std::tuple<Args...> m_subfields;
    std::vector<char> m_null_field_string;

    MultiFieldZ(RowZ *row, Args&&... subfields): m_row(row), m_subfields(subfields...) {
        init();
        m_null_field_string.assign(max_size(), 0);
    }

    /*
     * NdMultiFields must be bound to rows, and so all copies are assumed to be done as
     * part of a row copy. This is in contrast to the FieldBase class, which may be
     * legitimately copied without reference to a row.
     */
    MultiFieldZ(const MultiFieldZ& other) :
            m_row(other.m_row ? other.m_row->m_child : nullptr),
            m_subfields(other.m_subfields){
        init();
    }

    MultiFieldZ& operator=(const MultiFieldZ& other) {
        MPI_ASSERT(is_comparable(other), "Shapes are incompatible");
        struct fn_t {
            const MultiFieldZ& other;
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
            void operator()(FieldBaseZ &f) {
                MPI_REQUIRE(!f.m_row, "Field mustn't already be attached to a Row");
                f.m_row_offset = m_row->add_field(&f);
                f.m_row = m_row;
            }
        };
        fn_t fn {m_row};
        tuple_utils::for_each(m_subfields, fn);
    }

    bool operator==(const MultiFieldZ &other) const {
        struct fn_t {
            bool m_and = true;
            void operator()(const FieldBaseZ &f1, const FieldBaseZ &f2) { m_and &= f1==f2; }
        };
        fn_t fn;
        tuple_utils::for_each_pair(m_subfields, other.m_subfields, fn);
        return fn.m_and;
    }

    bool is_comparable(const MultiFieldZ &other) const {
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




#endif //M7_NDMULTIFIELDZ_H