//
// Created by Robert John Anderson on 2020-08-01.
//

#ifndef M7_CONCURRENTSTACK_H
#define M7_CONCURRENTSTACK_H


#include "ConcurrentLinkedList.h"

/**
 * First in, first out threadsafe stack implementation
 */
template<typename T>
struct ConcurrentStack : ConcurrentLinkedList<T> {

    using ConcurrentLinkedList<T>::m_next;

    bool pop(T &payload) {
        ConcurrentLinkedListNode<T> *first_node = nullptr;
#pragma omp atomic capture
        {
            first_node = m_next;
            m_next = first_node ? first_node->m_next : nullptr;
        }
        if (first_node) payload = m_next->m_payload;
        return first_node;
    }

};


#endif //M7_CONCURRENTSTACK_H
