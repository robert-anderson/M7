//
// Created by Robert J. Anderson on 1/25/22.
//

#ifndef M7_COMPOSITEFIELD_H
#define M7_COMPOSITEFIELD_H

#include <M7_lib/util/Tuple.h>

#include "M7_lib/hdf5/DatasetTransaction.h"
#include "FieldBase.h"

struct CompositeFieldBase {

    str_t prefix(str_t base, str_t prefix);

    void save(const hdf5::NodeWriter& nw, const str_t& name, uint_t max_nitem_per_op, bool this_rank) const {
        save_fn(nw, name, max_nitem_per_op, this_rank);
    }

    void save(const hdf5::NodeWriter& nw, const str_t& name, bool this_rank) const {
        save(nw, name, hdf5::c_default_max_nitem_per_op, this_rank);
    }

    void save(const hdf5::NodeWriter& nw, bool this_rank) const {
        // todo: get proper name
        save(nw, "", this_rank);
    }

    /*
     * functors for implementing the composite analogues of single-Field methods
     */
protected:

    virtual void save_fn(const hdf5::NodeWriter& nw, const str_t& name, uint_t max_nitem_per_op, bool this_rank) const = 0;

    struct EqFn {
        bool m_and = true;
        template<typename T>
        void operator()(const T &v1, const T &v2) {m_and &= (v1 == v2);}
    };


    struct IneqFn {
        bool m_gt_found_first = false;
        bool m_lt_found_first = false;
        template<typename T>
        void operator()(const T &f1, const T &f2) {
            if (f1 < f2) m_lt_found_first = !m_gt_found_first;
            else if (f1 > f2) m_gt_found_first = !m_lt_found_first;
        }
    };

    struct ToStringFn {
        str_t m_str;
        template<typename T>
        void operator()(const T &v) { m_str += v.to_string() + " "; }
    };

    struct ClearFn {
        template<typename T>
        void operator()(T &v) { v.clear();}
    };

    struct IsClearFn {
        bool m_and = true;
        template<typename T>
        void operator()(const T &v) { m_and &= v.is_clear();}
    };

    struct HashFn {
        hash::digest_t m_hash = 0;
        template<typename T>
        void operator()(const T &v) { m_hash ^= v.hash(); }
    };

    struct SaveFn {
        const hdf5::NodeWriter& m_nw;
        const str_t& m_name;
        const uint_t m_max_nitem_per_op;
        const bool m_this_rank;
        template<typename T>
        void operator()(const T &v) { v.save(m_nw, m_name, m_max_nitem_per_op, m_this_rank); }
    };

};


template<typename ...Args>
struct CompositeField : CompositeFieldBase {
    /**
     * can be made up of FieldBase subclasses or other CompositeFields
     */
    std::tuple<Args&...> m_refs;
    CompositeField(Args&... refs) : m_refs(refs...) {
        static_assert(sizeof...(Args)>0, "Composite field must have a non-zero number of components");
    }

    const Row* row() const {
        return get<0>().row();
    }

    bool belongs_to_row(const Row* row) const {
        return get<0>().belongs_to_row(row);
    }

    bool belongs_to_row() const {
        return belongs_to_row(row());
    }

    bool operator==(const CompositeField &other) const {
        EqFn fn;
        tuple::foreach(m_refs, other.m_refs, fn);
        return fn.m_and;
    }

    bool operator!=(const CompositeField &other) const {
        return !(*this == other);
    }

    bool operator<(const CompositeField &other) const{
        IneqFn fn;
        tuple::foreach(m_refs, other.m_refs, fn);
        return fn.m_lt_found_first;
    }

    bool operator<=(const CompositeField &other) const {
        IneqFn fn;
        tuple::foreach(m_refs, other.m_refs, fn);
        return fn.m_lt_found_first || (!fn.m_lt_found_first && !fn.m_gt_found_first);
    }

    bool operator>(const CompositeField &other) const{
        return !(*this<=other);
    }

    bool operator>=(const CompositeField &other) const {
        return !(*this<other);
    }

    str_t to_string() const {
        ToStringFn fn;
        tuple::foreach(m_refs, fn);
        return fn.m_str;
    }

    void clear() {
        ClearFn fn;
        tuple::foreach(m_refs, fn);
    }

    bool is_clear() const {
        IsClearFn fn;
        tuple::foreach(m_refs, fn);
        return fn.m_and;
    }

    hash::digest_t hash() const {
        HashFn fn;
        tuple::foreach(m_refs, fn);
        return fn.m_hash;
    }

    template<uint_t ifield>
    const typename std::tuple_element<ifield, std::tuple<Args...>>::type &get() const {
        return std::get<ifield>(m_refs);
    }

    template<uint_t ifield>
    typename std::tuple_element<ifield, std::tuple<Args...>>::type &get() {
        return std::get<ifield>(m_refs);
    }

protected:
    void save_fn(const hdf5::NodeWriter& nw, const str_t& name, uint_t max_nitem_per_op, bool this_rank) const override {
        SaveFn fn {nw, name, max_nitem_per_op, this_rank};
        tuple::foreach(m_refs, fn);
    }
};

template<typename ...Args>
static std::ostream &operator<<(std::ostream &os, const CompositeField<Args...> &v) {
    os << v.to_string();
    return os;
}


#endif //M7_COMPOSITEFIELD_H
