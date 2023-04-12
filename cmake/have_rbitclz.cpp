//
// Created by rja on 08/04/23.
//

#include <cstdint>

int main() {
    uint64_t n = 1;
    uint64_t count;
    asm("rbit %1, %1\nclz %0, %1": "=r" (count): "r" (n));
    return count;
}