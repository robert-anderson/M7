//
// Created by rja on 08/04/23.
//

#include <cstdint>

int main() {
    uint64_t n64 = 1;
    uint64_t res64;
    asm("tzcnt %1, %0;": "=r" (res64): "r" (n64));
    return n64;
}