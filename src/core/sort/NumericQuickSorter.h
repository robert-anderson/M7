//
// Created by Robert John Anderson on 2020-08-02.
//

#ifndef M7_NUMERICQUICKSORTER_H
#define M7_NUMERICQUICKSORTER_H

#include "QuickSorter.h"

template<typename T>
class NumericQuickSorter : public QuickSorter<T> {
public:
    NumericQuickSorter(const size_t &row_dsize, T *array) : QuickSorter<T>(row_dsize, array) {}

protected:
    int cmp(const T *row1, const T *row2) override {
        return *row1==*row2 ? 0 : (*row1<*row2?-1:1);
    }
};

template<typename T>
class NumericQuickSorter<std::complex<T>> : public QuickSorter<std::complex<T>> {
    NumericQuickSorter(const size_t &row_dsize, std::complex<T> *array) : QuickSorter<std::complex<T>>(row_dsize, array) {}
protected:
    int cmp(const T *row1, const T *row2) override {
        return *row1==*row2 ? 0 : (std::abs(*row1)<std::abs(*row2)?-1:1);
    }
};

#endif //M7_NUMERICQUICKSORTER_H
