//
// Created by Robert J. Anderson on 24/08/2021.
//

#include <M7_lib/basis/AbelianGroup.h>

#include "DecodedMbf.h"

void decoded_mbf::SimpleContainer::clear() {
    m_inds.clear();
}

bool decoded_mbf::SimpleContainer::empty() {
    return m_inds.empty();
}
