//
// Created by Robert John Anderson on 2020-02-09.
//

#include "BitfieldHasher.h"

#if 0
size_t BitfieldHasher::operator()(const BitfieldNew &key) const {
    size_t hash = 0ul;
    for (size_t idataword = 0ul; idataword < key.m_ndataword; ++idataword) {
        hash ^= std::hash<size_t>()(key.get_dataword(idataword));
    }
    return hash;
}
#endif
