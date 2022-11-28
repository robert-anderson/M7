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
    const str_t m_delimiters;
    str_t m_work_line;
public:

    CsvFileReader(const str_t& fname, str_t delimiters=", ()", uint_t iline=0ul);

    bool next(strv_t& tokens);
};

class NumericCsvFileReader : public CsvFileReader {
    static str_t c_allowed_chars;
    bool valid_numeric(const str_t& token);
    bool valid_numeric(const strv_t& tokens);

public:
    const uint_t m_ncolumn;
    NumericCsvFileReader(const str_t& fname, uint_t ncolumn,
                         str_t delimiters=", ()", uint_t iline=0ul);

    bool next(strv_t& tokens);

    typedef strv_t::const_iterator c_iter_token_t;

    static bool parsable_as(const str_t& str, uint_t&);
    static bool parsable_as(const str_t& str, int&);
    static bool parsable_as(const str_t& /*str*/, double&);
    static bool parsable_as(const str_t& /*str*/, float&);

    static void parse(const str_t& str, uint_t& v);
    static void parse(const str_t& str, long& v);
    static void parse(const str_t& str, int& v);
    static void parse(const str_t& str, double & v);
    static void parse(const str_t& str, float & v);

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
    static void parse(c_iter_token_t begin, c_iter_token_t end, v_t<T>& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        v.clear();
        for (auto it = begin; it!=end; ++it) {
            v.emplace_back();
            parse(it, it+1, v.back());
        }
    }

    template<typename T>
    static void parse(c_iter_token_t begin, c_iter_token_t end, v_t<std::complex<T>>& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        v.clear();
        for (auto it = begin; it != end; ++it) {
            v.emplace_back();
            parse(*it, arith::real_ref(v.back()));
            if (++it == end) return;
            parse(*it, arith::imag_ref(v.back()));
        }
    }

    static uint_t ncolumn(const str_t& fname, std::function<uint_t(const str_t&)> iline_fn);
};

#endif //M7_CSVFILEREADER_H
