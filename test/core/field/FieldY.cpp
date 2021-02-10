//
// Created by rja on 07/02/2021.
//

#include "gtest/gtest.h"

TEST(FieldY, Test){

#if 0
    std::vector<defs::data_t> m_buffer(100, 0);
    RowY row;
    row.m_dbegin = m_buffer.data();

    NdMultiFieldY<2, FermionOnvFieldY, BosonOnvFieldY> field(&row, {5, 4}, {2}, {8});

    std::cout << field.to_string() << std::endl;
    //field.get<1>().set(4);

    //std::cout << field.get<1>().to_string() << std::endl;

    std::string s = "hello";
    bit_utils::set(s[1], 3);
    std::cout << s << std::endl;
#endif

}