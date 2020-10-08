//
// Created by Robert John Anderson on 2020-08-02.
//

#include <src/core/sample/PRNG.h>
#include "gtest/gtest.h"
#include "src/core/sort/NumericQuickSorter.h"


TEST(QuickSorter, Test){
    std::vector<int> array = {3, 1, 4, 5, 7,34, 123, 523, 2, 1, 234, 34, 2};
    const size_t n = array.size();
    NumericQuickSorter<int> qsr(1, array.data());
    qsr.sort(n);
    for (size_t i=0ul; i<n; ++i) std::cout << array[i] << " ";
    std::cout << std::endl;

//    typedef size_t T;
//    struct alignas(defs::ncacheline_byte) Aligned {
//        T v;
//    };
//
//    std::vector<Aligned, AlignedAllocator<size_t, defs::ncacheline_byte>> v(n);
//    PRNG prng(123, n);
//    for (size_t i=0ul; i<n; ++i) v[i].v = prng.draw_uint()%100;
//
//    auto cmp = [&](const size_t& irow, const size_t& jrow){
//        ASSERT(irow<n && jrow<n)
//        const auto& iv = v[irow].v;
//        const auto& jv = v[jrow].v;
//        if (iv==jv) return 0;
//        if (iv<jv) return -1;
//        else return 1;
//    };
//
//    QuickSorter<Aligned> qs(1, v.data(), cmp);
//    for (size_t i=0ul; i<n; ++i) std::cout << v[i].v << " ";
//    std::cout << std::endl;
//    qs.sort(n);
//    for (size_t i=0ul; i<n; ++i) std::cout << v[i].v << " ";
//    std::cout << std::endl;
//
//
//    //auto t = omp_get_wtime();
//    //std::cout << nthread << " threads: " <<  omp_get_wtime()-t << std::endl;
}