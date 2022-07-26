//
// Created by rja on 26/07/22.
//

#include "version.h"

std::string get_version() {
    std::string str(VERSION);
    auto begin = str.cbegin();
    // strip off initial quote mark if present
    if (str.front()=='\'') ++begin;
    auto end = str.cend();
    // strip off final quote mark if present
    if (str.back()=='\'') --end;
    return "#"+std::string(begin, end);
}
