//
// Created by anderson on 1/25/22.
//

#ifndef M7_COMPOSITEFIELD_H
#define M7_COMPOSITEFIELD_H

#include "src/core/util/utils.h"

struct CompositeFieldBase {

    std::string prefix(std::string base, std::string prefix);

    /*
     * functors for implementing the composite analogues of single-Field methods
     */
protected:
    struct EqFn {
        bool m_and = true;
        template<typename T>
        void operator()(const T &v1, const T &v2) {m_and &= (v1 == v2);}
    };

    struct ToStringFn {
        std::string m_str;
        template<typename T>
        void operator()(const T &v) { m_str += v.to_string() + " "; }
    };
};


template<typename ...Args>
struct CompositeField : CompositeFieldBase {
    /**
     * can be made up of FieldBase descendants or other CompositeFields
     */
    std::tuple<Args&...> m_refs;
    CompositeField(Args&... refs) : m_refs(refs...) {}

    bool operator==(const CompositeField &other) const {
        EqFn fn;
        tuple_utils::for_each_pair(m_refs, other.m_refs, fn);
        return fn.m_and;
    }

    bool operator!=(const CompositeField &other) const {
        return !(*this == other);
    }

    std::string to_string() const {
        ToStringFn fn;
        tuple_utils::for_each(m_refs, fn);
        return fn.m_str;
    }

    template<size_t ifield>
    const typename std::tuple_element<ifield, std::tuple<Args...>>::type &get() const {
        return std::get<ifield>(m_refs);
    }

    template<size_t ifield>
    typename std::tuple_element<ifield, std::tuple<Args...>>::type &get() {
        return std::get<ifield>(m_refs);
    }
};


#endif //M7_COMPOSITEFIELD_H
