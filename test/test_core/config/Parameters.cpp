//
// Created by rja on 23/06/2021.
//

#include <M7_lib/conf/ConfComponents.h>
#include "gtest/gtest.h"

namespace parameters_test {

    struct Section1 : conf_components::Section {
        conf_components::Param<std::vector<size_t>> m_some_numbers;
        conf_components::Param<std::vector<size_t>> m_some_unspecified_numbers;
        conf_components::Param<std::string> m_a_string;

        struct SubSection1 : conf_components::Section {
            conf_components::Param<size_t> m_a_number;
            conf_components::Param<double> m_a_float;
            conf_components::Param<size_t> m_another_number;

            SubSection1(conf_components::Group *parent) :
                    conf_components::Section(parent, "subsection1", "minor options"),
                    m_a_number(this, "a_number", 0ul, "bla blah minor number"),
                    m_a_float(this, "a_float", 0.0, "bla blah minor float"),
                    m_another_number(this, "another_number", 6ul, "bla blah another minor number") {}
        };

        struct SubSection2 : conf_components::Section {
            conf_components::Param<size_t> m_a_number;

            SubSection2(conf_components::Group *parent) :
                    conf_components::Section(parent, "subsection2", "different minor options"),
                    m_a_number(this, "a_number", 0ul, "bla blah different minor number") {}
        };

        SubSection1 m_subsection1;
        SubSection2 m_subsection2;

        Section1(conf_components::Group *parent) :
                conf_components::Section(parent, "section1", "parameter nodes relating to section1"),
                m_some_numbers(this, "some_numbers", {3, 4, 6, 9}, "blah blah numbers"),
                m_some_unspecified_numbers(this, "some_unspecified_numbers", {3, 4, 6},
                                           "these numbers are not in the YAML file, and so the default value should be assigned"),
                m_a_string(this, "a_string", {}, "blah blah a string"),
                m_subsection1(this), m_subsection2(this) {}
    };

    struct Section2 : conf_components::Section {
        Section2(conf_components::Group *parent) : conf_components::Section(parent, "section2", "blah blah empty section") {}
    };

    struct Section3 : conf_components::Section {
        conf_components::Param <std::vector<double>> m_some_numbers;
        conf_components::Param<bool> m_some_flag;

        Section3(conf_components::Group *parent) :
                conf_components::Section(parent, "section3", "blah blah final section"),
                m_some_numbers(this, "some_numbers", {3, 4, 6}, "blah blah final section numbers"),
                m_some_flag(this, "some_flag", false, "blah blah final section flag") {}
    };

    struct TestDocument : conf_components::Document {
        Section1 m_section1;
        Section2 m_section2;
        Section3 m_section3;

        explicit TestDocument(const yaml::File *yf) :
                conf_components::Document(yf, "FCIQMC options", "options describing the behavior of something"),
                m_section1(this), m_section2(this), m_section3(this) {}
    };
}

TEST(Parameters, ParsingYaml) {
    yaml::File yf(defs::assets_root + "/yaml_test/example.yaml");
    ASSERT_TRUE(yf.exists("section1"));
    ASSERT_TRUE(yf.exists("section1.some_numbers"));
    ASSERT_EQ(yf.get_as<defs::inds>("section1.some_numbers"), defs::inds({2, 3, 5, 1}));
    ASSERT_EQ(yf.get_as<std::string>("section1.a_string"), "this is just a string");
    ASSERT_TRUE(yf.exists("section1.subsection1"));
    ASSERT_TRUE(yf.exists("section1.subsection1.a_number"));
    ASSERT_EQ(yf.get_as<size_t>("section1.subsection1.a_number"), 78);
    ASSERT_TRUE(yf.exists("section1.subsection1.a_float"));
    ASSERT_EQ(yf.get_as<double>("section1.subsection1.a_float"), 4.5);
    ASSERT_TRUE(yf.exists("section2"));
    ASSERT_TRUE(yf.exists("section3"));
    ASSERT_FALSE(yf.exists("section4"));
}

TEST(Parameters, DefaultValues) {
    using namespace parameters_test;
    TestDocument doc(nullptr);
    ASSERT_EQ(doc.m_section1.m_some_unspecified_numbers.get(), defs::inds({3, 4, 6}));
    auto str = doc.help_string();
    std::cout << str << std::endl;
}

TEST(Parameters, ParameterGroups) {
    yaml::File yf(defs::assets_root + "/yaml_test/example.yaml");
    using namespace parameters_test;
    TestDocument doc(&yf);

    auto invalid = doc.invalid_file_key();
    std::cout << "invalid key: " << invalid << std::endl;
    ASSERT_EQ(doc.m_section1.m_some_numbers.get(), defs::inds({2, 3, 5, 1}));
}