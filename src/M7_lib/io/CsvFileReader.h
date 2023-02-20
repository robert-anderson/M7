//
// Created by Robert J. Anderson on 12/3/21.
//

#ifndef M7_CSVFILEREADER_H
#define M7_CSVFILEREADER_H

#include <algorithm>

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/util/Parse.h>

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
    bool all_float_parseable(const strv_t& tokens);

public:
    const uint_t m_ncolumn;
    NumericCsvFileReader(const str_t& fname, uint_t ncolumn, str_t delimiters=", ()", uint_t iline=0ul);

    bool next(strv_t& tokens);

    typedef strv_t::const_iterator c_iter_token_t;

    template<typename T>
    static void parse(c_iter_token_t begin, c_iter_token_t /*end*/, T& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        auto success = parse::checked(*begin, v);
        DEBUG_ONLY(success);
        DEBUG_ASSERT_TRUE(success, logging::format("string \"{}\" could not be parsed", *begin));
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
