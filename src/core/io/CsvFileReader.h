//
// Created by anderson on 12/3/21.
//

#ifndef M7_CSVFILEREADER_H
#define M7_CSVFILEREADER_H

#include <algorithm>
#include <util/utils.h>
#include <parallel/MPIAssert.h>
#include "FileReader.h"

class CsvFileReader : public FileReader {
protected:
    const std::string m_delimiters;
    std::string m_work_line;
public:

    CsvFileReader(const std::string& fname, std::string delimiters=", ", size_t iline=0ul);

    bool next(std::vector<std::string>& tokens);
};

class NumericCsvFileReader : public CsvFileReader {
    static std::string c_allowed_chars;
    bool valid_numeric(const std::string& token);
    bool valid_numeric(const std::vector<std::string>& tokens);

public:
    const size_t m_ncolumn;
    NumericCsvFileReader(const std::string& fname, size_t ncolumn,
                         std::string delimiters=", ", size_t iline=0ul);

    bool next(std::vector<std::string>& tokens);

    typedef std::vector<std::string>::const_iterator c_iter_token_t;

    static bool parsable_as(const std::string& str, size_t& v);
    static bool parsable_as(const std::string& str, int& v);
    static bool parsable_as(const std::string& str, double& v);
    static bool parsable_as(const std::string& str, float& v);

    static void parse(const std::string& str, size_t& v);
    static void parse(const std::string& str, long& v);
    static void parse(const std::string& str, int& v);
    static void parse(const std::string& str, double & v);
    static void parse(const std::string& str, float & v);

    template<typename T>
    static void parse(c_iter_token_t begin, c_iter_token_t end, T& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        parse(*begin, v);
    }

    template<typename T>
    static void parse(c_iter_token_t begin, c_iter_token_t end, std::complex<T>& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        parse(*begin, consts::real_ref(v));
        if (begin+1==end) return;
        parse(*(begin+1), consts::imag_ref(v));
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
            parse(*it, consts::real_ref(v.back()));
            if (++it == end) return;
            parse(*it, consts::imag_ref(v.back()));
        }
    }

    static size_t ncolumn(const std::string& fname, std::function<size_t(const std::string&)> iline_fn);
};

#endif //M7_CSVFILEREADER_H
