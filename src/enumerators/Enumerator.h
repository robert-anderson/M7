//
// Created by Robert John Anderson on 2020-01-31.
//

#ifndef M7_ENUMERATOR_H
#define M7_ENUMERATOR_H

#include <vector>


template <typename result_T>
class Enumerator {
protected:
    Enumerator *m_subsequent = nullptr;
    Enumerator *m_current = this;
    Enumerator(){}
    Enumerator(Enumerator *subsequent){
        set_subsequent(subsequent);
    }

public:

    virtual inline bool next(result_T &result) {
        auto tmp = m_current->next_element(result);
        if (!tmp){
            if (!m_current->has_subsequent()) return false;
            else m_current = m_current->m_subsequent;
            return m_current->next(result);
        }
        return true;
    }
    virtual inline bool next(result_T &result, size_t &i) {
        ++i;
        return next(result);
    }

    void set_subsequent(Enumerator* subsequent){
        m_subsequent = subsequent;
    }
    bool has_subsequent(){
        return m_subsequent!= nullptr;
    }
    std::vector<result_T> enumerate() {
        std::vector<size_t> result{};
        size_t item;
        while (next(item)) {
            result.push_back(item);
        }
        return result;
    }
    size_t count(){
        size_t n = 0ul;
        result_T tmp;
        while(next(tmp)) n++;
        return n;
    }
    bool has_fewer_than_n_elements(size_t nmax) {
        size_t n = 0ul;
        result_T tmp;
        while(next(tmp)) {
            n++;
            if (n==nmax) return false;
        }
        return true;
    }

private:
    virtual inline bool next_element(result_T &result){}
};

#endif //M7_ENUMERATOR_H
