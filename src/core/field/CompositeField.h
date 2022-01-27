//
// Created by anderson on 1/25/22.
//

#ifndef M7_COMPOSITEFIELD_H
#define M7_COMPOSITEFIELD_H

#include <tuple>
#include "Fields.h"


struct CompositeFieldBase {

    /*
     * functors for implementing the composite analogues of single-Field methods
     */
protected:
    struct EqFn {
        bool m_and = true;
        template<typename T>
        void operator()(const T &v1, const T &v2) {m_and &= (v1 == v2);}
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
