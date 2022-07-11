//
// Created by Robert J. Anderson on 23/06/2021.
//

#include "gtest/gtest.h"
#include <M7_lib/conf/ConfComponents.h>

namespace parameters_test {

    using namespace conf_components;
    
    struct Section1 : Section {
        Param<v_t<uint_t>> m_some_numbers;
        Param<v_t<uint_t>> m_some_unspecified_numbers;
        Param<str_t> m_a_string;

        struct SubSection1 : Section {
            Param<uint_t> m_a_number;
            Param<double> m_a_float;
            Param<uint_t> m_another_number;

            SubSection1(conf_components::Group *parent) :
                    Section(parent, "subsection1", "minor options", conf_components::Implicit),
                    m_a_number(this, "a_number", 0ul, "bla blah minor number"),
                    m_a_float(this, "a_float", 0.0, "bla blah minor float"),
                    m_another_number(this, "another_number", 6ul, "bla blah another minor number") {}
        };

        struct SubSection2 : Section {
            Param<uint_t> m_a_number;

            SubSection2(conf_components::Group *parent) :
                    Section(parent, "subsection2", "different minor options"),
                    m_a_number(this, "a_number", 0ul, "bla blah different minor number") {}
        };

        SubSection1 m_subsection1;
        SubSection2 m_subsection2;

        Section1(conf_components::Group *parent) :
                Section(parent, "section1", "parameter nodes relating to section1"),
                m_some_numbers(this, "some_numbers", {3, 4, 6, 9}, "blah blah numbers"),
                m_some_unspecified_numbers(this, "some_unspecified_numbers", {3, 4, 6},
                                           "these numbers are not in the YAML file, and so the default value should be assigned"),
                m_a_string(this, "a_string", {}, "blah blah a string"),
                m_subsection1(this), m_subsection2(this) {}
    };

    struct Section2 : Section {
        Section2(conf_components::Group *parent) : Section(parent, "section2", "blah blah empty section") {}
    };

    struct Section3 : Section {
        Param <v_t<double>> m_some_numbers;
        Param<bool> m_some_flag;

        Section3(conf_components::Group *parent) :
                Section(parent, "section3", "blah blah final section"),
                m_some_numbers(this, "some_numbers", {3, 4, 6}, "blah blah final section numbers"),
                m_some_flag(this, "some_flag", false, "blah blah final section flag") {}
    };

    struct TestDocument : conf_components::Document {
        Section1 m_section1;
        Section2 m_section2;
        Section3 m_section3;

        explicit TestDocument(const str_t& fname) :
                conf_components::Document(fname, "options describing the behavior of something"),
                m_section1(this), m_section2(this), m_section3(this) {}
    };
}

#if 0
TEST(Parameters, ParsingYaml) {
    conf_components::Document yf(PROJECT_ROOT"/assets/yaml_test/example.yaml", "");
    ASSERT_TRUE(yf.exists("section1"));
    ASSERT_TRUE(yf.exists("section1.some_numbers"));
    ASSERT_EQ(yf.get_as<uintv_t>("section1.some_numbers"), uintv_t({2, 3, 5, 1}));
    ASSERT_EQ(yf.get_as<str_t>("section1.a_string"), "this is just a string");
    ASSERT_TRUE(yf.exists("section1.subsection1"));
    ASSERT_TRUE(yf.exists("section1.subsection1.a_number"));
    ASSERT_EQ(yf.get_as<uint_t>("section1.subsection1.a_number"), 78);
    ASSERT_TRUE(yf.exists("section1.subsection1.a_float"));
    ASSERT_EQ(yf.get_as<double>("section1.subsection1.a_float"), 4.5);
    ASSERT_TRUE(yf.exists("section2"));
    ASSERT_TRUE(yf.exists("section3"));
    ASSERT_TRUE(yf.exists("section4"));
    ASSERT_FALSE(yf.exists("section5"));
}

TEST(Parameters, DefaultValues) {
    using namespace parameters_test;
    TestDocument doc(nullptr);
    ASSERT_EQ(doc.m_section1.m_some_unspecified_numbers.get(), uintv_t({3, 4, 6}));
    auto str = doc.help_string();
}

TEST(Parameters, ParameterGroups) {
    conf_components::Document yf(PROJECT_ROOT"/assets/yaml_test/example.yaml");
    using namespace parameters_test;
    TestDocument doc(&yf);
    auto invalid = doc.invalid_file_key();
    ASSERT_EQ(doc.m_section1.m_some_numbers.get(), uintv_t({2, 3, 5, 1}));
}


#endif


TEST(Parameters, New){
    using namespace parameters_test;
    TestDocument f(PROJECT_ROOT"/assets/yaml_test/example.yaml");
    f.validate();
    f.print_help();
    f.log();

//    auto i = Node(n, "section2").nchild_in_file();
//    std::cout << i << std::endl;
}