//
// Created by rja on 09/02/2021.
//

#ifndef M7_MULTIFIELD_H
#define M7_MULTIFIELD_H

#include "Row.h"
#include "FieldBase.h"
#include "src/core/nd/NdSelector.h"


template<typename ...Args>
struct MultiField {
    Row *m_row = nullptr;
    std::tuple<Args...> m_subfields;
    std::vector<char> m_null_field_string;

    MultiField(Row *row, Args&&... subfields): m_subfields(subfields...) {
        add_to_row(row);
        m_null_field_string.assign(max_size(), 0);
    }

    /*
     * NdMultiFields must be bound to rows, and so all copies are assumed to be done as
     * part of a row copy. This is in contrast to the FieldBase class, which may be
     * legitimately copied without reference to a row.
     */
    MultiField(const MultiField& other) :
            m_row(other.m_row ? other.m_row->m_child : nullptr),
            m_subfields(other.m_subfields){
        add_to_row(m_row);
    }

    MultiField& operator=(const MultiField& other) {
        MPI_ASSERT(is_comparable(other), "Shapes are incompatible");
        struct fn_t {
            const MultiField& other;
            void operator()(FieldBase &f1, const FieldBase &f2) { f1=f2; }
        };
        fn_t fn;
        tuple_utils::for_each_pair(m_subfields, other.m_subfields, fn);
        return *this;
    }

    void add_to_row(Row* row) {
        if (!row) return;
        MPI_REQUIRE(!m_row, "MultiField must not be already associated with a row");
        struct fn_t {
            Row* m_row;
            void operator()(FieldBase &f) {
                f.add_to_row(m_row);
            }
        };
        m_row = row;
        fn_t fn {m_row};
        tuple_utils::for_each(m_subfields, fn);
    }

    bool operator==(const MultiField &other) const {
        struct fn_t {
            bool m_and = true;
            void operator()(const FieldBase &f1, const FieldBase &f2) { m_and &= f1 == f2; }
        };
        fn_t fn;
        tuple_utils::for_each_pair(m_subfields, other.m_subfields, fn);
        return fn.m_and;
    }

    bool operator!=(const MultiField &other) const {
        return !(*this==other);
    }

    bool operator<(const MultiField &other) const {
        struct fn_t {
            bool m_gt_found_first = false;
            bool m_lt_found_first = false;
            void operator()(const FieldBase &f1, const FieldBase &f2) {
                if (f1 < f2) m_lt_found_first = m_gt_found_first ? false : true;
                else if (f1 > f2) m_gt_found_first = m_lt_found_first ? false : true;
            }
        };
        fn_t fn;
        tuple_utils::for_each_pair(m_subfields, other.m_subfields, fn);
        return fn.m_lt_found_first;
    }

    bool operator<=(const MultiField &other) const {
        struct fn_t {
            bool m_gt_found_first = false;
            bool m_lt_found_first = false;
            void operator()(const FieldBase &f1, const FieldBase &f2) {
                if (f1 < f2) m_lt_found_first = m_gt_found_first ? false : true;
                else if (f1 > f2) m_gt_found_first = m_lt_found_first ? false : true;
            }
        };
        fn_t fn;
        tuple_utils::for_each_pair(m_subfields, other.m_subfields, fn);
        return fn.m_lt_found_first || (!fn.m_lt_found_first && !fn.m_gt_found_first);
    }

    bool operator>(const MultiField &other) const {
        return !(*this<=other);
    }

    bool operator>=(const MultiField &other) const {
        return !(*this<other);
    }

    bool is_comparable(const MultiField &other) const {
        struct fn_t {
            bool m_and = true;
            void operator()(const FieldBase &f1, const FieldBase &f2) { m_and &= f1.is_comparable(f2); }
        };
        fn_t fn;
        tuple_utils::for_each_pair(m_subfields, other.m_subfields, fn);
        return fn.m_and;
    }

    bool max_size() const {
        struct fn_t {
            size_t m_max = 0ul;
            void operator()(const FieldBase &f) { if (f.m_size > m_max) m_max = f.m_size; }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
        return fn.m_max;
    }

    void zero() {
        struct fn_t {
            void operator()(FieldBase &f) { f.zero(); }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
    }

    bool is_zero() const {
        struct fn_t {
            bool m_and = true;
            void operator()(const FieldBase &f) { m_and &= f.is_zero(); }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
        return fn.m_and;
    }

    bool is_added_to_row() const {
        struct fn_t {
            bool m_and = true;
            void operator()(const FieldBase &f) { m_and &= f.belongs_to_row(); }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
        return fn.m_and;
    }

    defs::hash_t hash() const {
        struct fn_t {
            defs::hash_t m_hash = 0;
            void operator()(const FieldBase &f) { m_hash ^= f.hash(); }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
        return fn.m_hash;
    }

public:
    std::string to_string() const {
        struct fn_t {
            std::string m_str;

            void operator()(const FieldBase &f) { m_str += f.to_string() + " "; }
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

template<typename ...Args>
static std::ostream &operator<<(std::ostream &os, const MultiField<Args...> &v) {
    os << v.to_string();
    return os;
}

#endif //M7_MULTIFIELD_H