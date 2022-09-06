/*
 * this file was generated programmatically and should not be modified by hand.
 * 
 * only a subset of spin signature pairs are valid for RDM processing, so the
 * array in this file maps the creation and annihilation spin signature to a spin
 * signature "case_id" for subsequent processing
 */
#ifndef M7_SPINSIGMAP_H
#define M7_SPINSIGMAP_H
#include "M7_lib/defs.h"
constexpr uint_t spinsig_map[16][16] = 
      {{0, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul},
       {~0ul, 1, 3, ~0ul, 7, ~0ul, ~0ul, ~0ul, 21, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul},
       {~0ul, 2, 4, ~0ul, 9, ~0ul, ~0ul, ~0ul, 23, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul},
       {~0ul, ~0ul, ~0ul, 5, ~0ul, 12, 15, ~0ul, ~0ul, 28, 35, ~0ul, 47, ~0ul, ~0ul, ~0ul},
       {~0ul, 6, 8, ~0ul, 10, ~0ul, ~0ul, ~0ul, 25, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul},
       {~0ul, ~0ul, ~0ul, 11, ~0ul, 13, 17, ~0ul, ~0ul, 30, 37, ~0ul, 49, ~0ul, ~0ul, ~0ul},
       {~0ul, ~0ul, ~0ul, 14, ~0ul, 16, 18, ~0ul, ~0ul, 32, 39, ~0ul, 51, ~0ul, ~0ul, ~0ul},
       {~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, 19, ~0ul, ~0ul, ~0ul, 44, ~0ul, 58, 63, ~0ul},
       {~0ul, 20, 22, ~0ul, 24, ~0ul, ~0ul, ~0ul, 26, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul},
       {~0ul, ~0ul, ~0ul, 27, ~0ul, 29, 31, ~0ul, ~0ul, 33, 41, ~0ul, 53, ~0ul, ~0ul, ~0ul},
       {~0ul, ~0ul, ~0ul, 34, ~0ul, 36, 38, ~0ul, ~0ul, 40, 42, ~0ul, 55, ~0ul, ~0ul, ~0ul},
       {~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, 43, ~0ul, ~0ul, ~0ul, 45, ~0ul, 60, 65, ~0ul},
       {~0ul, ~0ul, ~0ul, 46, ~0ul, 48, 50, ~0ul, ~0ul, 52, 54, ~0ul, 56, ~0ul, ~0ul, ~0ul},
       {~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, 57, ~0ul, ~0ul, ~0ul, 59, ~0ul, 61, 67, ~0ul},
       {~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, 62, ~0ul, ~0ul, ~0ul, 64, ~0ul, 66, 68, ~0ul},
       {~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, ~0ul, 69}};

#endif //M7_SPINSIGMAP_H
