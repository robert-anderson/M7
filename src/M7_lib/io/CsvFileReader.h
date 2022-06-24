//
// Created by Robert J. Anderson on 12/3/21.
//

#ifndef M7_CSVFILEREADER_H
#define M7_CSVFILEREADER_H

#include <algorithm>

#include <M7_lib/parallel/MPIAssert.h>

#include "FileReader.h"

class CsvFileReader : public FileReader {
protected:
    const std::string m_delimiters;
    std::string m_work_line;
public:

    CsvFileReader(const std::string& fname, std::string delimiters=", ", uint_t iline=0ul);

    bool next(std::vector<std::string>& tokens);
};

class NumericCsvFileReader : public CsvFileReader {
    static std::string c_allowed_chars;
    bool valid_numeric(const std::string& token);
    bool valid_numeric(const std::vector<std::string>& tokens);

public:
    const uint_t m_ncolumn;
    NumericCsvFileReader(const std::string& fname, uint_t ncolumn,
                         std::string delimiters=", ", uint_t iline=0ul);

    bool next(std::vector<std::string>& tokens);

    typedef std::vector<std::string>::const_iterator c_iter_token_t;

    static bool parsable_as(const std::string& str, uint_t&);
    static bool parsable_as(const std::string& str, int&);
    static bool parsable_as(const std::string& /*str*/, double&);
    static bool parsable_as(const std::string& /*str*/, float&);

    static void parse(const std::string& str, uint_t& v);
    static void parse(const std::string& str, long& v);
    static void parse(const std::string& str, int& v);
    static void parse(const std::string& str, double & v);
    static void parse(const std::string& str, float & v);

    template<typename T>
    static void parse(c_iter_token_t begin, c_iter_token_t /*end*/, T& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        parse(*begin, v);
    }

    template<typename T>
    static void parse(c_iter_token_t begin, c_iter_token_t end, std::complex<T>& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        parse(*begin, arith::real_ref(v));
        if (begin+1==end) return;
        parse(*(begin+1), arith::imag_ref(v));
    }

    template<typename T>
    static void parse(c_iter_token_t begin, c_iter_token_t end, std::vector<T>& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        v.clear();
        for (auto it = begin; it!=end; ++it) {
            v.emplace_back();
            parse(it, it+1, v.back());
        }
    }

    template<typename T>
    static void parse(c_iter_token_t begin, c_iter_token_t end, std::vector<std::complex<T>>& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        v.clear();
        for (auto it = begin; it != end; ++it) {
            v.emplace_back();
            parse(*it, arith::real_ref(v.back()));
            if (++it == end) return;
            parse(*it, arith::imag_ref(v.back()));
        }
    }

    static uint_t ncolumn(const std::string& fname, std::function<uint_t(const std::string&)> iline_fn);
};

#endif //M7_CSVFILEREADER_H
