//
// Created by rja on 21/10/2020.
//

#include "gtest/gtest.h"
#include "src/core/data/Table.h"
#include "src/core/data/BitsetField.h"
#include "src/core/data/DeterminantField.h"
#include "src/core/data/BosonOnvField.h"


TEST(Table_ND, Test){

    TableX t;
    defs::inds data(100, 0);
    t.m_data = (char*)data.data();
    //BitsetFieldX<0> field(&t, {}, 10, "asdasdsadas");
    //DeterminantFieldX<0> dets(&t, {}, 5, "asdasdsadas");
    BosonOnvField<0> perms(&t, {}, 8, "zddsaas");

    auto view = perms(0);
    std::cout << view.to_string() << std::endl;

//    std::vector<size_t> vec(6);
//    auto p = (size_t*)((char*)vec.data()+1);
//    p[0] = ~0ul;
//    for (auto i: vec) std::cout << i << std::endl;


//    Table t;
//    std::vector<double> vec(50, 9.99);
//    t.m_data = (char*)vec.data();
//    NumericArrayField<double, 1, 2> matrices(&t, {6}, {2, 3}, "some matrices");
//    matrices(0, 0)(0, 0) = 0;
//    matrices(0, 0)(0, 1) = 1;
//    matrices(0, 0)(0, 2) = 2;
//    matrices(0, 0)(1, 0) = 3;
//    matrices(0, 0)(1, 1) = 4;
//    matrices(0, 0)(1, 2) = 5;
//    matrices(0, 2) = matrices(0, 0);
//    matrices.raw_view(0, 0).second;

    //for (auto item: vec) std::cout << item << std::endl;
}
