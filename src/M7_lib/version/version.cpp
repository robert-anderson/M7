//
// Created by rja on 26/07/22.
//

#include "version.h"

std::string get_version() {
    std::string str(VERSION);
    return "#"+str;
}

std::string get_compilation_timestamp() {
    return {COMPILATION_TIMESTAMP};
}
