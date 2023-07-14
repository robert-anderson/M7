//
// Created by anderson on 14/07/2023.
//

#include <queue>
#include <M7_lib/util/Integer.h>
#include <M7_lib/util/Bit.h>
#include <M7_lib/util/Hash.h>
#include "test_core/defs.h"

class IndexSet {
    std::map<uint_t, uint_t> m_map;

public:

    IndexSet(const uintv_t& v) {
        for (auto& i: v) insert(i);
    }

    void insert(uint_t i) {
        const auto iword = i / (sizeof(uint_t) * CHAR_BIT);
        const auto ibit = i % (sizeof(uint_t) * CHAR_BIT);
        auto it = m_map.find(iword);
        if (it == m_map.cend()) {
            uint_t content = 0ul;
            bit::set(content, ibit);
            m_map.insert({iword, content});
        }
        else bit::set(it->second, ibit);
    }

    bool operator[](uint_t i) const {
        const auto iword = i / (sizeof(uint_t) * CHAR_BIT);
        auto it = m_map.find(iword);
        if (it == m_map.cend()) return false;
        const auto ibit = i % (sizeof(uint_t) * CHAR_BIT);
        return bit::get(it->second, ibit);
    }

    void clear() {
        m_map.clear();
    }

    bool empty() {
        return m_map.empty();
    }

    static void intersection(const IndexSet& a, const IndexSet& b, IndexSet& c) {
        c.clear();
        if (&a == &c || &b == &c) return;
        auto it_a = a.m_map.cbegin();
        auto it_b = b.m_map.cbegin();
        while (it_a != a.m_map.cend() && it_b != b.m_map.cend()) {
            if (it_a->first < it_b->first) ++it_a;
            else if (it_a->first > it_b->first) ++it_b;
            else {
                // same group of inds
                const auto isect = it_a->second & it_b->second;
                // only add isect to result if it's not completely clear
                if (isect) c.m_map.insert({it_a->first, isect});
                ++it_a; ++it_b;
            }
        }
    }

    IndexSet& intersect_with(const IndexSet& other) {
        auto it = m_map.begin();
        while (it != m_map.end()) {
            auto oit = other.m_map.find(it->first);
            if (oit != other.m_map.cend()) {
                it->second &= oit->second;
                if (it->second) {
                    ++it;
                    continue;
                }
            }
            it = m_map.erase(it);
        }
        return *this;
    }

    template<typename fn_t>
    void foreach(const fn_t& fn) const {
        functor::assert_prototype<void(uint_t)>(fn);
        for (auto& pair: m_map) {
            auto work = pair.second;
            while (work) {
                const auto ibit = bit::next_setbit(work);
                fn(pair.first * (sizeof(uint_t) * CHAR_BIT) + ibit);
            }
        }
    }

    std::set<uint_t> to_std_set() const {
        std::set<uint_t> s;
        auto fn = [&s](uint_t i) {s.insert(i);};
        foreach(fn);
        return s;
    }

};

TEST(IndexSet, ToStdSet) {
    const uint_t nelem = 10000;
    auto v = hash::unique_in_range(0, nelem, 0, 3*nelem);
    const std::set<uint_t> s(v.cbegin(), v.cend());
    IndexSet is(v);
    const auto s_from_is = is.to_std_set();
    ASSERT_EQ(s, s_from_is);
}