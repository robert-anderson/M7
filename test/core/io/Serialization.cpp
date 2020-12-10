//
// Created by rja on 09/12/2020.
//

#include "gtest/gtest.h"
#include "src/core/io/Serialization.h"

namespace serialization_test {
    struct TestSerializable : Serializable {
        double d = 12113425.4;
        float f = 1231.4;
        std::vector<char> v = {1,12,3,43,53,6,65,7,8};
        TestSerializable(){
            store(d);
            store(f);
            store(v);
        }
    };
}

TEST(Serialization, Test){
    serialization_test::TestSerializable ts;
    ASSERT_EQ(ts.dsize(), 6);
    std::cout << ts.format_hash()<< std::endl;

    Serializer s(100);
    s.set_output("tmp.bin");
    s.set_input("tmp.bin");
    ts.save_write(s);
    ts.v[4] = 22;
    for (auto c: ts.v) {std::cout << (int)c << " ";} std::cout << std::endl;
    ts.read_load(s);
    for (auto c: ts.v) {std::cout << (int)c << " ";} std::cout << std::endl;

}