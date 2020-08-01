//
// Created by Robert John Anderson on 2020-08-01.
//

#ifndef M7_CONCURRENTLINKEDLIST_H
#define M7_CONCURRENTLINKEDLIST_H

#include <cstddef>
#include <iostream>

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
        auto tmp = node->m_next;
        if (!tmp) return;
        node->m_next = tmp->m_next;
        if (tmp==m_last) m_last = node;
        delete tmp;
    }

    bool is_empty(){
        return m_next;
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
        ConcurrentLinkedListNode<T> *old_last_node;
        auto new_node = new CarrierNode<T>(payload);
#pragma omp atomic capture
        {
            old_last_node = m_last;
            m_last = new_node;
        }
        old_last_node->m_next = new_node;
        return old_last_node;
    }

    void print() {
        for (CarrierNode<size_t>* node = this->m_next; node; node = node->m_next)
            std::cout << node->m_payload << " ";
        std::cout << std::endl;
    }

    /*
     * for debugging only! not a constant-time operation
     */
    size_t size() const {
        size_t n = 0ul;
        for (CarrierNode<size_t>* node = this->m_next; node; node = node->m_next) ++n;
        return n;
    }

};


#endif //M7_CONCURRENTLINKEDLIST_H
