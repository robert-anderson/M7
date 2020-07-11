//
// Created by rja on 10/07/2020.
//

#include <src/core/io/FileIterator.h>
#include <src/core/util/defs.h>
#include "gtest/gtest.h"

TEST(FileIterator, Test) {
    FileIterator iterator(defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP");
    std::string line;
    while (iterator.next(line)){
        std::cout << line << std::endl;
    }
}