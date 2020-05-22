//
// Created by Robert John Anderson on 2020-03-29.
//

#ifndef M7_HASHTABLE_H
#define M7_HASHTABLE_H


#include <exception>
#include <stack>
#include "src/core/util/defs.h"
#include "HashMap.h"

#if 0
class KeyError : public std::exception {
    virtual const char *what() const throw() {
        return "Key does not exist in hash table";
    }
};

/*
 * A HashMap that stores its own key-value pairs
 */

template<typename key_T, typename value_T>
class HashTable : public HashMap {
    std::vector<std::pair<key_T, value_T>> m_pairs;
    std::stack<size_t> m_free;
public:

    HashTable(const size_t &nbucket) : HashMap(nbucket), m_pairs(nbucket) {
        for (size_t i = 0ul; i < nbucket; ++i) {
            m_free.push(i);
        }
    }

    key_T get_key(const size_t &key_index) const override {
        return m_pairs[key_index].first;
    }

    void set_key(const size_t &key_index, const key_T &key) override {
        m_pairs[key_index].first = key;
    }

    HashTable(HashTable &old, const size_t &nbucket) :
        HashMap(old, nbucket), m_pairs(nbucket) {
        m_pairs.assign(old.m_pairs.begin(), old.m_pairs.end());
        m_free = old.m_free;
    }

    const value_T get(const key_T &key) {
        auto ipair = HashMap<key_T>::lookup(key);
        if (ipair != ~0ul) return m_pairs[ipair].second;
        throw KeyError();
    }

    void set(const key_T &key, const value_T &value) {
        auto ipair = HashMap<key_T>::lookup(key);
        if (ipair != ~0ul) m_pairs[ipair].second = value;
        else {
            size_t ipair = m_free.top();
            m_free.pop();
            HashMap<key_T>::insert(key, ipair);
            m_pairs[ipair].second = value;
        }
    }

    size_t remove(const key_T &key) override {
        auto ipair = HashMap<key_T>::remove(key);
        m_free.push(ipair);
        return ipair;
    }
};



#endif //M7_HASHTABLE_H
#endif //M7_HASHTABLE_H
