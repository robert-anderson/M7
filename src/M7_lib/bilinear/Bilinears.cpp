//
// Created by Robert J. Anderson on 11/08/2021.
//

#include "Bilinears.h"
#include "M7_lib/connection/OpSig.h"

OpSig bilinears::parse_exsig(const str_t& string) {
    REQUIRE_TRUE_ALL(string.size() == 1 || string.size() == 4, "invalid exsig string specification");
    if (string.size() == 1) {
        uint_t rank = string::parse_decimal_digit(string.c_str());
        REQUIRE_LE_ALL(rank, opsig::c_nop_mask_frm, "number of fermion operators exceeds limit");
        return opsig::frm(rank);
    }
    uint_t nfrm_cre = string::parse_decimal_digit(string.c_str());
    REQUIRE_LE_ALL(nfrm_cre, opsig::c_nop_mask_frm, "number of fermion creation operators exceeds limit");
    uint_t nfrm_ann = string::parse_decimal_digit(string.c_str()+1);
    REQUIRE_LE_ALL(nfrm_ann, opsig::c_nop_mask_frm, "number of fermion annihilation operators exceeds limit");
    uint_t nbos_cre = string::parse_decimal_digit(string.c_str()+2);
    REQUIRE_LE_ALL(nbos_cre, opsig::c_nop_mask_bos, "number of boson creation operators exceeds limit");
    uint_t nbos_ann = string::parse_decimal_digit(string.c_str()+3);
    REQUIRE_LE_ALL(nbos_ann, opsig::c_nop_mask_bos, "number of boson annihilation operators exceeds limit");
    return {{nfrm_cre, nfrm_ann}, {nbos_cre, nbos_ann}};
}

v_t<OpSig> bilinears::parse_exsigs(const strv_t& strings) {
    v_t<OpSig> out;
    for (auto &string: strings) out.push_back(parse_exsig(string));
    DEBUG_ASSERT_EQ(out.size(), strings.size(),
                    "output should have the same number of exsigs as specification");
    return out;
}
