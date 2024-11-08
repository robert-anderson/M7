//
// Created by rja on 08/04/23.
//

#include <cstdint>

int main() {
    uint32_t n32 = 0;
    uint64_t n64 = 0;
    uint32_t res32;
    uint64_t res64;
    res32 = __builtin_popcount(n32);
    res64 = __builtin_popcountl(n64);
    return n32 || n64;
}