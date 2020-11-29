//
// Created by rja on 29/11/2020.
//

#ifndef M7_QUICKSORTER_H
#define M7_QUICKSORTER_H


#include <src/defs.h>
#include <functional>

struct Quicksorter {

    defs::inds m_inds;
    typedef std::function<bool(const size_t &, const size_t &)> comp_t;
    comp_t m_comp_fn;

    Quicksorter(comp_t comp_fn);

    void sort(const size_t &hwm);

    bool is_sorted(const size_t &hwm);

private:
    void swap(size_t ii1, size_t ii2);

    size_t partition(size_t iilo, size_t iihi);

    void qs(size_t iilo, size_t iihi);

};


#endif //M7_QUICKSORTER_H
