//
// Created by rja on 24/08/2021.
//

#include "DecodedMbf.h"
#include "basis/AbelianGroup.h"

void decoded_mbf::SimpleContainer::clear() {
    m_inds.clear();
}

bool decoded_mbf::SimpleContainer::empty() {
    return m_inds.empty();
}