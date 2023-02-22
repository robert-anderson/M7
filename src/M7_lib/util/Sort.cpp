//
// Created by anderson on 22/02/2023.
//

#include <cstring>
#include <M7_lib/parallel/MPIAssert.h>
#include "Sort.h"

void sort::reorder(void *data, uint_t element_size, const uintv_t &order) {
    const auto size = element_size * order.size();
    v_t<buf_t> copy(size);
    std::memcpy(copy.data(), data, size);
    auto dst = reinterpret_cast<buf_t*>(data);
    for (auto i: order){
        DEBUG_ASSERT_LT(i, order.size(), "index OOB");
        std::memcpy(dst, copy.data() + element_size * i, element_size);
        dst += element_size;
    }
}
