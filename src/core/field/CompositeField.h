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
    CompositeField(Args &&... children) : m_children(std::forward<Args>(children)...) {
    }
};


#endif //M7_COMPOSITEFIELD_H
