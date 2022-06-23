//
// Created by rja on 15/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/Tuple.h"

namespace tuple_test {
    struct PlusOneFn {
        template<typename T>
        void operator()(T &v) { v+=1;}
        void operator()(std::string &v) {v.back()++;};
    };
    struct PlusOneChkFn {
        template<typename T>
        void operator()(const T &v1, const T &v2) {
            ASSERT_EQ(v1 + 1, v2);
        }
        void operator()(const std::string &v1, const std::string &v2) {
            auto v1_cpy = v1;
            v1_cpy.back()+=1;
            ASSERT_EQ(v1_cpy, v2);
        };
    };

}

TEST(Tuple, ForeachSingle) {
    std::tuple<std::string, int, float> tup = {"example", 4, -1.5};
    std::tuple<std::string, int, float> tup_chk = {"examplf", 5, -0.5};
    tuple_test::PlusOneFn fn;
    tuple::foreach(tup, fn);
    ASSERT_EQ(tup, tup_chk);
}

TEST(Tuple, ForeachPair) {
    std::tuple<std::string, int, float> tup = {"example", 4, -1.5};
    std::tuple<std::string, int, float> tup_chk = {"examplf", 5, -0.5};
    tuple_test::PlusOneChkFn fn;
    tuple::foreach(tup, tup_chk, fn);
}