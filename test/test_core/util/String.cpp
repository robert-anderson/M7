//
// Created by rja on 12/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/String.h"


TEST(UtilString, SplitLine) {
    const std::string line = "0.5000000000 1 1 2 2";
    auto tokens = utils::string::split(line, ' ');
}

TEST(UtilString, JoinAndSplit) {
    const std::string line = " this is   an   example   string   ";
    auto tokens = utils::string::split(line, ' ');
    ASSERT_EQ(tokens.size(), 5);
    auto joinder = utils::string::join(tokens, " ");
    // splitting will eliminate consecutive occurrences of the delimiter
    ASSERT_EQ("this is an example string", joinder);
}

TEST(UtilString, Tokenize) {
    const std::string line = " this is   an,   example   string   ";
    auto tokens = utils::string::split(line, " ,");
    ASSERT_EQ(tokens.size(), 5);
}