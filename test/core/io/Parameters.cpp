//
// Created by rja on 23/06/2021.
//

#include "gtest/gtest.h"
#include "src/core/io/Parameters.h"

namespace parameters_test {


    struct Section1 : config::Section {
        config::Param<std::vector<size_t>> m_some_numbers;
        config::Param<std::string> m_a_string;
        struct SubSection1 : config::Section {
            config::Param<size_t> m_a_number;
            config::Param<double> m_a_float;
            config::Param<size_t> m_another_number;
            SubSection1(config::Group *parent):
                    config::Section(parent, "subheading1", "minor options"),
                    m_a_number(this, "a_number", "bla blah minor number"),
                    m_a_float(this, "a_float", "bla blah minor float"),
                    m_another_number(this, "another_number", "bla blah another minor number"){}
        };
        struct SubSection2 : config::Section {
            config::Param<size_t> m_a_number;
            SubSection2(config::Group *parent):
                    config::Section(parent, "subheading2", "different minor options"),
                    m_a_number(this, "a_number", "bla blah different minor number"){}
        };
        SubSection1 m_subsection1;
        SubSection2 m_subsection2;

        Section1(config::Group *parent) :
                config::Section(parent, "heading1", "parameter nodes relating to heading1"),
                m_some_numbers(this, "some_numbers", "blah blah numbers", {3, 4, 6, 9}, config::check::size_eq<>()),
                m_a_string(this, "a_string", "blah blah a string"),
                m_subsection1(this), m_subsection2(this){}
    };

    struct Section2: config::Section {
        Section2(config::Group *parent) : config::Section(parent, "heading2", "blah blah empty section"){}
    };

    struct Section3: config::Section {
        config::Param<std::vector<double>> m_some_numbers;
        Section3(config::Group *parent) :
            config::Section(parent, "heading3", "blah blah final section"),
            m_some_numbers(this, "some_numbers", "blah blah final section numbers", {3, 4, 6}){}
    };

    struct TestDocument : config::Document {
        Section1 m_section1;
        Section2 m_section2;
        Section3 m_section3;

        TestDocument(const yaml::File *yf) :
                config::Document(yf, "parameter document", "options describing the behavior of something"),
                m_section1(this), m_section2(this), m_section3(this) {}
    };
}

TEST(Parameters, ParsingYaml) {
    yaml::File yf(defs::assets_root + "/yaml_test/example.yaml");
    ASSERT_TRUE(yf.exists("heading1"));
    ASSERT_TRUE(yf.exists("heading1.some_numbers"));
    ASSERT_EQ(yf.get_as<defs::inds>("heading1.some_numbers"), defs::inds({2, 3, 5, 1}));
    ASSERT_EQ(yf.get_as<std::string>("heading1.a_string"), "this is just a string");
    ASSERT_TRUE(yf.exists("heading1.subheading1"));
    ASSERT_TRUE(yf.exists("heading1.subheading1.a_number"));
    ASSERT_EQ(yf.get_as<size_t>("heading1.subheading1.a_number"), 78);
    ASSERT_TRUE(yf.exists("heading1.subheading1.a_float"));
    ASSERT_EQ(yf.get_as<double>("heading1.subheading1.a_float"), 4.5);
    ASSERT_TRUE(yf.exists("heading2"));
    ASSERT_TRUE(yf.exists("heading3"));
    ASSERT_FALSE(yf.exists("heading4"));
}

TEST(Parameters, DefaultValues) {
    using namespace parameters_test;
    TestDocument doc(nullptr);
    ASSERT_EQ(doc.m_section1.m_some_numbers.get(), defs::inds({3, 4, 6}));
}

TEST(Parameters, ParameterGroups) {
    yaml::File yf(defs::assets_root + "/yaml_test/example.yaml");
    using namespace parameters_test;
    TestDocument doc(nullptr);
    //TestDocument doc(&yf);

    auto invalid = doc.invalid_file_key();
    std::cout << "invalid key: " << invalid << std::endl;
    ASSERT_EQ(doc.m_section1.m_some_numbers.get(), defs::inds({2, 3, 5, 1}));
}