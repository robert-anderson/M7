//
// Created by Robert John Anderson on 2020-08-01.
//

#ifndef M7_CONCURRENTLINKEDLIST_H
#define M7_CONCURRENTLINKEDLIST_H

#include <cstddef>
#include <iostream>
#include "src/core/util/defs.h"

template<typename T>
struct CarrierNode;

template<typename T>
struct ConcurrentLinkedListNode {
    CarrierNode<T> *m_next = nullptr;
};

template<typename T>
struct CarrierNode : public ConcurrentLinkedListNode<T> {
    T m_payload;
    explicit CarrierNode(const T& payload): m_payload(payload){}
};


template<typename T>
struct ConcurrentLinkedList : ConcurrentLinkedListNode<T> {
    using ConcurrentLinkedListNode<T>::m_next;
    ConcurrentLinkedListNode<T> *m_last = this;

    void delete_after(ConcurrentLinkedListNode<T> *node) {
        ASSERT(node)
        auto tmp = node->m_next;
        if (!tmp) {
            throw std::runtime_error("can't delete the next node since it does not exist");
        }
        node->m_next = tmp->m_next;
        if (tmp==m_last) m_last = node;
        delete tmp;
    }

    bool is_empty() const {
        return !m_next;
    }

    void clear() {
        while (m_next) delete_after(this);
    }

    ~ConcurrentLinkedList() {
        clear();
    }

    /*
     * returns the node pointing to the appended node
     */
    ConcurrentLinkedListNode<T> *append(const T &payload) {
        ASSERT((m_last==this) || !is_empty())
        ConcurrentLinkedListNode<T> *old_last_node;
        auto new_node = new CarrierNode<T>(payload);
#pragma omp atomic capture
        {
            old_last_node = m_last;
            m_last = new_node;
        }
        ASSERT(old_last_node)
        old_last_node->m_next = new_node;
        ASSERT(!is_empty())
        ASSERT(m_last!=this)
        return old_last_node;
    }

    std::vector<T> to_vector() const {
        std::vector<T> out;
        for (CarrierNode<T>* node = this->m_next; node; node = node->m_next)
            out.push_back(node->m_payload);
        return out;
    }

    /*
     * for debugging only! not a constant-time operation
     */
    size_t size() const {
        size_t n = 0ul;
        for (CarrierNode<T>* node = this->m_next; node; node = node->m_next) ++n;
        return n;
    }

};


#endif //M7_CONCURRENTLINKEDLIST_H
