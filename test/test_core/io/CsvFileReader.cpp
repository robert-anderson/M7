//
// Created by anderson on 12/3/21.
//

#include "src/core/io/CsvFileReader.h"
#include "gtest/gtest.h"

TEST(CsvFileReader, Fcidump){
    CsvFileReader file_reader(defs::assets_root+"/RHF_Cr2_12o12e/FCIDUMP");
    std::vector<std::string> split_line;
    file_reader.next(split_line);
    std::vector<std::string> chk_line;
    chk_line = {"&FCI", "NORB=", "12", "NELEC=12", "MS2=0"};
    ASSERT_EQ(split_line,  chk_line);
    for (size_t iline=1ul; iline<123; ++iline) {
        file_reader.next(split_line);
    }
    chk_line = {"-0.001368384767972724", "4", "3", "9", "9"};
    ASSERT_EQ(split_line,  chk_line);
}

TEST(NumericCsvFileReader, Fcidump){
    NumericCsvFileReader file_reader(defs::assets_root+"/RHF_Cr2_12o12e/FCIDUMP", 5);
    std::vector<std::string> split_line;
    file_reader.next(split_line);
    std::vector<std::string> chk_line;
    chk_line = {"0.5192717990625102", "1", "1", "1", "1"};
    ASSERT_EQ(split_line,  chk_line);
    for (size_t iline=5ul; iline<123; ++iline) {
        file_reader.next(split_line);
    }
    chk_line = {"-0.001368384767972724", "4", "3", "9", "9"};
    ASSERT_EQ(split_line,  chk_line);
}

TEST(NumericCsvFileReader, ParseScientific){
    double v;
    NumericCsvFileReader::parse("1.e6", v);
    ASSERT_FLOAT_EQ(v, 1000000.0);
    NumericCsvFileReader::parse("-0.4e3", v);
    ASSERT_FLOAT_EQ(v, -400.0);
    NumericCsvFileReader::parse("+3141.596e-3", v);
    ASSERT_FLOAT_EQ(v, 3.141596);
    NumericCsvFileReader::parse("+0.03141596E+2", v);
    ASSERT_FLOAT_EQ(v, 3.141596);
}

TEST(NumericCsvFileReader, ParseSingleReal){
    const size_t nreal = 1;
    std::vector<std::string> line = {"-0.001368384767972724", "4", "3", "9", "9"};
    defs::inds inds;
    defs::inds chk_inds = {4, 3, 9, 9};
    NumericCsvFileReader::parse(line.cbegin()+nreal, line.cend(), inds);
    ASSERT_EQ(inds, chk_inds);
    std::vector<double> real_values;
    std::vector<double> chk_real_values = {-0.001368384767972724};
    NumericCsvFileReader::parse(line.cbegin(), line.cbegin()+nreal, real_values);
    ASSERT_EQ(real_values, chk_real_values);
    std::vector<std::complex<double>> complex_values;
    std::vector<std::complex<double>> chk_complex_values = {{-0.001368384767972724, 0.0}};
    NumericCsvFileReader::parse(line.cbegin(), line.cbegin()+nreal, complex_values);
    ASSERT_EQ(complex_values, chk_complex_values);
}

TEST(NumericCsvFileReader, ParseMultiReal){
    const size_t nreal = 5;
    std::vector<std::string> line = {"-0.00134", "12.22", "-1.33", "4.5463", "-9.012", "4", "3", "9", "9"};
    defs::inds inds;
    defs::inds chk_inds = {4, 3, 9, 9};
    NumericCsvFileReader::parse(line.cbegin()+nreal, line.cend(), inds);
    ASSERT_EQ(inds, chk_inds);
    std::vector<double> real_values;
    std::vector<double> chk_real_values = {-0.00134, 12.22, -1.33, 4.5463, -9.012};
    NumericCsvFileReader::parse(line.cbegin(), line.cbegin()+nreal, real_values);
    ASSERT_EQ(real_values, chk_real_values);
    std::vector<std::complex<double>> complex_values;
    std::vector<std::complex<double>> chk_complex_values = {{-0.00134, 12.22}, {-1.33, 4.5463}, {-9.012, 0.0}};
    NumericCsvFileReader::parse(line.cbegin(), line.cbegin()+nreal, complex_values);
    ASSERT_EQ(complex_values, chk_complex_values);
}