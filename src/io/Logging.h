//
// Created by Robert John Anderson on 2020-01-04.
//

#include <iostream>
#include "src/parallel/MPIWrapper.h"

namespace logger {

    enum level{debug, info, warning};

    void write(std::string msg, size_t indent=0, level l=info){
        if (mpi::i_am_root()) {
            for (size_t i=0ul; i<indent; ++i) std::cout << "  ";
            std::cout << msg << std::endl;
        }
    }
}