//
// Created by anderson on 1/25/22.
//

#ifndef M7_COMPOSITEFIELD_H
#define M7_COMPOSITEFIELD_H

#include <tuple>

template<typename ...Args>
struct CompositeField {
    /**
     * can be made up of FieldBase descendants or other CompositeFields
     */
    std::tuple<Args...> m_children;
    CompositeField(Args &&... children) : m_children(std::forward<Args>(children)...) {}





    template<size_t ifield>
    const typename std::tuple_element<ifield, std::tuple<Args...>>::type &get() const {
        return std::get<ifield>(m_children);
    }

    template<size_t ifield>
    typename std::tuple_element<ifield, std::tuple<Args...>>::type &get() {
        return std::get<ifield>(m_children);
    }
};


#endif //M7_COMPOSITEFIELD_H
