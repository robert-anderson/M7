//
// Created by rja on 19/08/2021.
//

#ifndef M7_EBDUMPFILEREADER_H
#define M7_EBDUMPFILEREADER_H

#include "HamiltonianFileReader.h"

struct EbdumpFileReader : HamiltonianFileReader {
    EbdumpFileReader(const std::string &fname): HamiltonianFileReader(fname, 3, false){
        REQUIRE_FALSE_ALL(m_spin_resolved, "spin resolved electron-boson dumps are not currently supported");
    }

};


#endif //M7_EBDUMPFILEREADER_H
