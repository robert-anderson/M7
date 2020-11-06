//
// Created by Robert John Anderson on 2020-01-31.
//

#ifndef M7_ENUMERATOR_H
#define M7_ENUMERATOR_H

#include <vector>

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

#endif //M7_ENUMERATOR_H
