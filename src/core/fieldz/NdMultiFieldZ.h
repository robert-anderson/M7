//
// Created by rja on 09/02/2021.
//

#ifndef M7_NDMULTIFIELDZ_H
#define M7_NDMULTIFIELDZ_H

#include "RowZ.h"
#include "src/core/nd/NdSequence.h"

struct NdFieldBaseZ {
    RowZ *m_row;
    NdFieldBaseZ(RowZ* row): m_row(row){}
};

template<size_t nind, typename ...Args>
struct NdMultiFieldZ : NdFieldBaseZ {
    NdSequence<nind> m_sequence;
    mutable typename NdSequence<nind>::Cursor m_inds;
    std::tuple<Args...> m_subfields;

private:
    bool oob() const {
        return static_cast<const FieldBaseZ &>(get<0>()).oob();
    }
public:

    NdMultiFieldZ(std::array<size_t, nind> shape, Args&&... subfields) :
            NdFieldBaseZ(nullptr), m_sequence(shape), m_inds(m_sequence.cursor()), m_subfields(std::move(subfields)...) {
    }

    NdMultiFieldZ(RowZ *row, std::array<size_t, nind> shape, Args&&... subfields) :
            NdFieldBaseZ(row), m_sequence(shape), m_inds(m_sequence.cursor()), m_subfields(std::move(subfields)...) {
        init();
    }

    /*
     * NdMultiFields must be bound to rows, and so all copies are assumed to be done as
     * part of a row copy. This is in contrast to the FieldBase class, which may be
     * legitimately copied without reference to a row.
     */
    NdMultiFieldZ(const NdMultiFieldZ& other) :
            NdFieldBaseZ(other.m_row->m_child),
            m_sequence(other.m_sequence),
            m_inds(other.m_inds),
            m_subfields(other.m_subfields){
        init();
    }

    NdMultiFieldZ& operator=(const NdMultiFieldZ& other) {
        MPI_ASSERT(other.m_sequence.m_format==m_sequence.m_format, "Shapes are incompatible");
        m_row = other.m_row;
        m_subfields = other.m_subfields;
        return *this;
    }

    void init() {
        if (!m_row) return;
        struct fn_t {
            RowZ* m_row;
            size_t m_nelement;
            void operator()(FieldBaseZ &f) {
                f.m_nelement = m_nelement;
                f.m_size = f.m_nelement * f.m_element_size;
                f.m_row_offset = m_row->add_field(&f);
                f.m_row = m_row;
                f.m_max_view_offset = f.m_row_offset + f.m_size;
            }
        };
        fn_t fn {m_row, m_sequence.m_format.nelement()};
        tuple_utils::for_each(m_subfields, fn);
        restart();
    }

    bool equals(const NdMultiFieldZ &other) const {
        struct fn_t {
            bool m_and = true;

            void operator()(const FieldBaseZ &f1, const FieldBaseZ &f2) { m_and &= f1.equals(f2); }
        };
        fn_t fn;
        tuple_utils::for_each_pair(m_subfields, other.m_subfields, fn);
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

    bool restart() const {
        m_inds.to_front();
        struct fn_t {
            void operator()(const FieldBaseZ &f) { f.m_view_offset = f.m_row_offset; }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
        return get<0>().oob();
    }

    bool step() const {
        m_inds++;
        struct fn_t {
            void operator()(const FieldBaseZ &f) { f.m_view_offset += f.m_element_size; }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
        if (get<0>().m_view_offset >= get<0>().m_max_view_offset) {
            restart();
            return false;
        }
        return true;
    }

    template<typename ...Inds>
    bool select(Inds... inds) {
        m_inds.to(inds...);
        struct fn_t {
            size_t m_ielement;

            void operator()(const FieldBaseZ &f) { f.m_view_offset = f.m_row_offset + f.m_element_size * m_ielement; }
        };
        fn_t fn;
        tuple_utils::for_each(m_subfields, fn);
        if (get<0>().m_view_offset >= get<0>().m_max_view_offset) {
            restart();
            return false;
        }
        return true;
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


template<size_t nind, typename fieldbase_t>
struct NdFieldZ :  NdMultiFieldZ<nind, fieldbase_t>{
    NdFieldZ(RowZ *row, std::array<size_t, nind> shape, fieldbase_t&& subfield):
            NdMultiFieldZ<nind, fieldbase_t>(row, shape, std::move(subfield)){}

    const fieldbase_t &operator()() const {
        return NdMultiFieldZ<nind, fieldbase_t>::template get<0>();
    }

    fieldbase_t &operator()() {
        return NdMultiFieldZ<nind, fieldbase_t>::template get<0>();
    }
};

#endif //M7_NDMULTIFIELDZ_H