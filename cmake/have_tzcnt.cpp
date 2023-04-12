//
// Created by rja on 08/04/23.
//

#include <cstdint>

int main() {
    uint32_t n32 = 0;
    uint64_t n64 = 0;
    uint32_t res32;
    uint64_t res64;
    asm("tzcnt %1, %0;": "=r" (res32): "r" (n32));
    asm("tzcntq %1, %0;": "=r" (res64): "r" (n64));
    return n32 && n64;
}
