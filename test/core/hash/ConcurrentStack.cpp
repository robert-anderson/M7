//
// Created by Robert John Anderson on 2020-08-03.
//

#include "gtest/gtest.h"
#include "src/core/hash/ConcurrentStack.h"

TEST(ConcurrentStack, Test){
    ConcurrentStack<size_t> stack;
    size_t i;
    ASSERT_FALSE(stack.pop(i));
    stack.append(4);
    stack.append(5);
    stack.append(8);
    ASSERT_TRUE(stack.pop(i));
    ASSERT_EQ(i, 4);
    ASSERT_TRUE(stack.pop(i));
    ASSERT_EQ(i, 5);
    ASSERT_TRUE(stack.pop(i));
    ASSERT_EQ(i, 8);
    ASSERT_FALSE(stack.pop(i));
    ASSERT_FALSE(stack.pop(i));
    ASSERT_TRUE(stack.is_empty());
    stack.append(4);
    ASSERT_FALSE(stack.is_empty());
    stack.clear();
    ASSERT_TRUE(stack.is_empty());
    stack.append(4);
    ASSERT_FALSE(stack.is_empty());
}