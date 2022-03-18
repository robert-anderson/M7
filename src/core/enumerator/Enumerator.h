//
// Created by Robert John Anderson on 2020-01-31.
//

#ifndef M7_ENUMERATOR_H
#define M7_ENUMERATOR_H

#include <vector>
#include <algorithm>
#include <defs.h>
#include <util/utils.h>
#include <parallel/MPIAssert.h>

template <typename result_t>
class Enumerator {
    typedef result_t enum_result_t;
protected:
    Enumerator *m_subsequent = nullptr;
    Enumerator *m_current = this;
    Enumerator(){}
    Enumerator(Enumerator *subsequent){
        set_subsequent(subsequent);
    }

public:

    virtual ~Enumerator(){}

    virtual bool next(result_t &result) {
        auto tmp = m_current->next_element(result);
        if (!tmp){
            if (!m_current->has_subsequent()) return false;
            else m_current = m_current->m_subsequent;
            return m_current->next(result);
        }
        return true;
    }

    virtual bool next(result_t &result, size_t &i) {
        ++i;
        return next(result);
    }

    void set_subsequent(Enumerator* subsequent){
        m_subsequent = subsequent;
    }
    bool has_subsequent(){
        return m_subsequent!= nullptr;
    }

    result_t default_result(){
        return result_t{};
    }

    std::vector<result_t> enumerate() {
        std::vector<result_t> results{};
        auto result = default_result();
        while (next(result)) {
            results.push_back(result);
        }
        return results;
    }

    size_t count(result_t& v){
        size_t n = 0ul;
        result_t tmp = v;
        while(next(tmp)) n++;
        return n;
    }

    bool has_fewer_than_n_elements(size_t nmax) {
        size_t n = 0ul;
        result_t tmp;
        while(next(tmp)) {
            n++;
            if (n==nmax) return false;
        }
        return true;
    }

private:
    virtual bool next_element(result_t &result) = 0;
};


namespace enums {
    struct Enumerator {
        const size_t m_n;
        const size_t m_r;
        const size_t m_nv;
        defs::inds m_v;

        Enumerator(size_t n, size_t r, size_t nv);

        virtual bool next() = 0;

        const size_t &operator[](const size_t &i) const;
    };

    static constexpr size_t c_cache_limit = 1ul<<20;
    /*
     * Enumerators are supposed to be as efficient as possible, but all overhead involved in
     * generation of their terms can be eliminated at the point of use by caching terms, this
     * is the role of this class.
     */
    template<typename enum_T>
    struct Cache {
        const size_t m_n, m_r;
        const defs::inds m_v;
        const size_t m_nv;
        size_t m_iv = ~0ul;
    private:
        std::vector<size_t> make() const {
            enum_T e(m_n, m_r);
            /*
             * check whether this is a situation in which the programmer should be encouraged to use
             * the enum_T to compute each term on the fly. This class is only to be used in situations
             * where the number of terms is small, and the number of iterations is expected to be large
             */
            if (e.m_nv>c_cache_limit)
                log::warn("Attempting to create a Cache of over {} terms", c_cache_limit);
            defs::inds out;
            while (e.next()) out.insert(out.end(), e.m_v.cbegin(), e.m_v.cend());
            return out;
        }

    public:
        Cache(size_t n, size_t r) : m_n(n), m_r(r), m_v(make()), m_nv(m_v.size() / m_r) {
            std::cout << m_nv << "  " << enum_T(m_n, m_r).m_nv << std::endl;
            REQUIRE_EQ(m_nv, enum_T(m_n, m_r).m_nv, "Number of terms generated does not meet expectation");
        }

        bool next() {
            ++m_iv;
            return m_iv < m_nv;
        }

        const size_t &operator[](const size_t &i) const {
            ASSERT(i < m_r);
            return m_v[m_iv * m_r + i];
        }
    };

    struct CombinationsWithRepetition : Enumerator {
        CombinationsWithRepetition(size_t n, size_t r);

        bool next() override;
    };

    struct CombinationsDistinct : Enumerator {
        std::string m_starting_bitmask, m_bitmask;
        bool m_allfound = false;
        CombinationsDistinct(size_t n, size_t r);

        bool next() override;
    };

    struct PermutationsWithRepetition : Enumerator {
        const defs::inds m_shape;

        static defs::inds make_shape(size_t n, size_t r) {
            defs::inds shape;
            shape.reserve(r);
            shape.assign(r, n);
            return shape;
        }
        PermutationsWithRepetition(size_t n, size_t r);
        PermutationsWithRepetition(const defs::inds& shape);

        bool next() override;
    };
}


#endif //M7_ENUMERATOR_H
