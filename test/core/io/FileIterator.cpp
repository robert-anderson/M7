//
// Created by rja on 10/07/2020.
//

#include <src/core/io/FileIterator.h>
#include <src/core/util/defs.h>
#include "gtest/gtest.h"


template<typename T, size_t nind>
struct SparseArrayFileIterator {//: public FileIterator {

    enum Arithmetic {
        None,
        Real,
        Complex
    };
    Arithmetic m_arithmetic = None;
    bool m_brackets = false;
    bool m_indsfirst = false;

    SparseArrayFileIterator(){}//FileIterator(){}

    /*
    bool set_kind(const std::string& line){
        ASSERT(m_arithmetic==None);
        size_t iopen = line.find('(');
        size_t iclose = line.find(')');
        if (iopen!=~0ul && iclose==~0ul)
            return false; // brackets are invalid
        if (iopen==~0ul && iclose!=~0ul)
            return false; // brackets are invalid

        std::cout << iopen-~0ul << std::endl;
        m_brackets = iopen;
    }*/

};

double read_double(char*& ptr){
    auto is_numeric = [](const char& c){
        return '0'<=c && c<='9';
    };
    auto is_partial_standard_float = [&](const char& c){
        return is_numeric(c) || c=='.' || c=='-';
    };
    auto is_partial_scientific = [&](const char& c){
        return is_partial_standard_float(c) || c=='e' || c=='E'|| c=='d'|| c=='D'|| c=='+';
    };
    auto is_divider = [](const char& c){
        return c==' ' || c==',' || c==')';
    };
    char* begin = nullptr;
    //for(; *ptr!=0; ptr++){
    while(*ptr!=0){
        if (!begin){
            if (is_partial_standard_float(*ptr)) begin = ptr;
            DBVAR(*ptr)
        }
        else {
            if (is_divider(*ptr)) {
                DBVAR(*ptr)
                return std::strtod(begin, &ptr);
            }
            else if (!is_partial_scientific(*ptr)) {
                begin = nullptr;
                DBVAR(*ptr)
            }
        }
        ptr++;
    }
    DBVAR((size_t)ptr)
    if (begin && begin!=ptr) {
        DBVAR(*ptr)
        return std::strtod(begin, &ptr);
    }
    else {
        DBVAR(*ptr)
        return std::numeric_limits<double>::max();
    }
}

char read(char*& p){
    char c = *p;
    for(; *p!=0; p++){}
    return c;
}

TEST(FileIterator, Test) {
    std::string line = "TREL=.TRUE.";
    char* p = line.begin().base();
    while (*p!=0){

        //DBVAR((size_t)p)
        //std::cout << read_double(p) << std::endl;
        std::cout << read(p) << std::endl;
        //DBVAR((size_t)p)


    }

//    FileIterator iterator(defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP");
//    std::string line;
//    while (iterator.next(line)){
//        std::cout << line << std::endl;
//        ASSERT(*line.end().base()==0)
//        char* p = line.begin().base();
//        while (*p!=0){
//            std::cout << p << "     " << read_double(p) << std::endl;
//        }
//    }
}