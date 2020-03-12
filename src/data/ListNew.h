//
// Created by rja on 12/03/2020.
//

#ifndef M7_LISTNEW_H
#define M7_LISTNEW_H

#include <src/defs.h>
#include "TableNew.h"

struct ListNew : public TableNew {
    defs::inds m_high_water_mark;

    explicit ListNew(size_t nsegment=1): TableNew(nsegment), m_high_water_mark(nsegment, 0){}

    virtual size_t push(const size_t &isegment=0) {
        size_t tmp;
#pragma omp atomic capture
        tmp = m_high_water_mark[isegment]++;
        if (tmp>=m_nrow_per_segment) throw std::runtime_error("Reached capacity of List");
        return tmp;
    }

    size_t push(const size_t &isegment, const size_t &nrow) {
        size_t tmp;
#pragma omp atomic capture
        tmp = m_high_water_mark[isegment] += nrow;
        if (tmp>=m_nrow_per_segment) throw std::runtime_error("Reached capacity of List");
        return tmp;
    }

    void zero() override{
        // TODO: no actual need to memset zero here, only included initially for clarity in debugging
        TableNew::zero();
        m_high_water_mark.assign(m_nsegment, 0);
    }

    //void print(size_t irank) const override;

};


#endif //M7_LISTNEW_H
