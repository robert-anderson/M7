//
// Created by anderson on 10/02/2022.
//

#include "gtest/gtest.h"
#include "src/core/table/BufferedFields.h"

TEST(BitsetField, SetFromInds) {

    size_t ibegin = 8;
    size_t iend = 9;

    typedef short T;
    T buf = 0;
    buf = ~buf;
    const size_t nbit_word = sizeof(T)*CHAR_BIT;
    buf = (buf >> T(nbit_word-iend)) & (buf << ibegin);

    for (size_t ibit=0; ibit<nbit_word; ++ibit)
        std::cout << (bit_utils::get(buf, ibit) ? '1':'0');
    std::cout << std::endl;

    //const size_t nbit = 100;
    // deliberately using non-standard container type
    //buffered::Bitset<short> field(nbit);
    //std::cout << field << std::endl;
    //field.set_range(0, 5);
}
