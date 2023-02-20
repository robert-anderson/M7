//
// Created by anderson on 20/02/2023.
//

#ifndef M7_PARSE_H
#define M7_PARSE_H

#include "M7_lib/defs.h"

namespace parse {
    static void with_except(const str_t& str, uint_t& v) {
        v = std::stoul(str);
    }

    static void with_except(const str_t& str, long& v) {
        v = std::stol(str);
    }

    static void with_except(const str_t& str, int& v) {
        v = std::stoi(str);
    }

    static void with_except(const str_t& str, double& v) {
        v = std::stod(str);
    }

    template<typename T>
    bool checked(const str_t& str, T& v) {
        try {
            with_except(str, v);
        }
        catch (const std::invalid_argument& ex){
            return false;
        }
        catch (const std::out_of_range& ex){
            return false;
        }
        return true;
    }
};


#endif //M7_PARSE_H
