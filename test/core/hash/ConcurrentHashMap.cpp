//
// Created by Robert John Anderson on 2020-07-30.
//

#include "gtest/gtest.h"
#include <omp.h>

template<typename T>
struct LinkedList {
    struct Node {
        T m_content;
        Node *m_next = nullptr;
    };
    Node *m_first = nullptr;
    Node *m_last = nullptr;

    void delete_after(const Node *node) {
        Node *tmp;
        if (!node) {
            if (m_first) {
                tmp = m_first->m_next;
                if (m_last == m_first) m_last = nullptr;
                delete m_first;
                m_first = tmp;
            }
        } else {
            if (node != m_last) {
                tmp = node->m_next;
                const_cast<Node *&>(node)->m_next = tmp->m_next;
                delete tmp;
            }
        }
    }

    void clear() {
        while (m_first) delete_after(nullptr);
    }

    ~LinkedList() {
        clear();
    }

    Node *append(const T &index) {
        Node *old_last_node = nullptr;
        auto new_node = new Node;
        new_node->m_content = index;
#pragma omp atomic capture
        {
            old_last_node = m_last;
            m_last = new_node;
        }
        if (old_last_node) old_last_node->m_next = new_node;
        else m_first = new_node;
        return new_node;
    }

    void print() {
        Node *node = m_first;
        while (node) {
            std::cout << node->m_content << " ";
            node = node->m_next;
        }
        std::cout << std::endl;
    }

};


struct HashMapBucket : public LinkedList<size_t>{
    LinkedList<LinkedList<size_t>::Node *> m_tombstone_prevs;

    void mark_for_delete_after(LinkedList<size_t>::Node *node){
        m_tombstone_prevs.append(node);
        if (node) node->m_next->m_content = ~0ul;
        else m_first->m_content = ~0ul;
    }

    void clear_tombstones() {
        auto tombstone_prev = m_tombstone_prevs.m_first;
        while (tombstone_prev) {
            delete_after(tombstone_prev->m_content);
            tombstone_prev = tombstone_prev->m_next;
        }
    }
};

template<typename T>
struct HashMap {
    std::vector<HashMapBucket> m_buckets;
    void clear_tombstones() {
        // pragma omp parallel for
        for (auto bucket : m_buckets) bucket.clear_tombstones();
    }
};

TEST(ConcurrentHashMap, Test) {
    HashMapBucket bucket;
    bucket.append(4);
    bucket.append(41);
    bucket.append(42);
    bucket.append(43);
    bucket.mark_for_delete_after(nullptr);
    bucket.print();
    bucket.clear_tombstones();
    bucket.print();
}



