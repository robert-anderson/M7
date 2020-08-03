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
    using ConcurrentLinkedList<T>::m_last;

    bool pop(T &payload) {
        CarrierNode<T> *first_node = nullptr;
//#pragma omp atomic capture
        {
            first_node = m_next;
            m_next = first_node ? first_node->m_next : nullptr;
        }
        if (first_node) {
            payload = first_node->m_payload;
            delete first_node;
            return true;
        } else {
            m_last=this;
            ASSERT(!m_next)
            return false;
        }
    }

};


#endif //M7_CONCURRENTSTACK_H
