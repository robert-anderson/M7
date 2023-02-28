//
// Created by rja on 12/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/String.h"


TEST(UtilString, SplitLine) {
    const str_t line = "0.5000000000 1 1 2 2";
    auto tokens = string::split(line, ' ');
}

TEST(UtilString, JoinAndSplit) {
    const str_t line = " this is   an   example   string   ";
    auto tokens = string::split(line, ' ');
    ASSERT_EQ(tokens.size(), 5);
    auto joinder = string::join(tokens, " ");
    // splitting will eliminate consecutive occurrences of the delimiter
    ASSERT_EQ("this is an example string", joinder);
}

TEST(UtilString, Tokenize) {
    const str_t line = " this is   an,   example   string   ";
    auto tokens = string::split(line, " ,");
    ASSERT_EQ(tokens.size(), 5);
}

TEST(UtilString, ToUpper) {
    str_t str = "tHis iS an eXample strIng 123";
    string::to_upper(str);
    ASSERT_EQ(str, "THIS IS AN EXAMPLE STRING 123");
}

TEST(UtilString, ToLower) {
    str_t str = "tHis iS an eXample strIng 123";
    string::to_lower(str);
    ASSERT_EQ(str, "this is an example string 123");
}