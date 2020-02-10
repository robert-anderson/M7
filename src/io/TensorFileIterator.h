//
// Created by Robert John Anderson on 2020-01-17.
//

#ifndef M7_TENSORFILEITERATOR_H
#define M7_TENSORFILEITERATOR_H

#include "FileIterator.h"
#include "../defs.h"
#include "../utils.h"
#include <iostream>
#include <assert.h>
#include <complex>

const std::string float_regex_string = R"([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)";
const std::string complex_regex_string = R"(\()"+float_regex_string+R"(\,)"+float_regex_string+R"(\))";
const std::string space_delimited_uint_regex_string = R"((^|\s)\d+($|\s))";

const std::regex float_regex(float_regex_string);
const std::regex complex_regex(complex_regex_string);
const std::regex inds_first_regex(space_delimited_uint_regex_string+".*"+float_regex_string);
const std::regex float_inds_regex(
        space_delimited_uint_regex_string+".*"+float_regex_string+"|"+
        float_regex_string+".*"+space_delimited_uint_regex_string);
const std::regex space_delimited_uint_regex(space_delimited_uint_regex_string);


template <typename T>
class TensorFileIterator : public FileIterator{
    static_assert(std::__is_static_castable<T, float>::value || std::__is_static_castable<T, std::complex<float>>::value);
public:
    static constexpr size_t m_nreal = std::__is_static_castable<T, float>::value?1:2;
    const size_t m_nind;
    const bool m_indsfirst;
private:
    const std::string m_data_regex_string;
    const std::regex m_data_regex;
public:

    const size_t m_nreal_given;
    const defs::inds m_shape;

    TensorFileIterator(const std::string &filename, const size_t &nind, const bool &indsfirst):
            FileIterator(filename, std::regex(data_regex_string(nind, indsfirst))),
            m_nind(nind), m_indsfirst(indsfirst), m_data_regex_string(data_regex_string(nind, indsfirst)),
            m_data_regex(std::regex(m_data_regex_string)), m_nreal_given(nreal_given(filename, m_data_regex)),
            m_shape(shape(filename, m_nind, m_nreal_given, m_indsfirst, m_data_regex)){}

    static std::string uint_space_list_regex_string(const size_t &n) {
        return R"((^|\s))"+string_utils::join(R"(\d+)", n, R"(\s+)")+R"(($|\s))";
    }

    static std::string data_regex_string(const size_t &nind, const bool &indsfirst){
        if (indsfirst) return uint_space_list_regex_string(nind)+R"(.*)"+float_regex_string;
        else return float_regex_string+R"(.*)"+uint_space_list_regex_string(nind);
    }

    bool next(defs::inds &inds, T &value) {
        std::string line;
        FileIterator::next(line);
        if (line.empty()) return 0;
        extract_line(line, m_nind, m_nreal_given, m_indsfirst, inds, value);
        return !line.empty();
    }

    static size_t nind(const std::string &filename){
        FileIterator fi(filename, float_inds_regex);
        std::string line{fi.next()};
        auto out{0ul};
        auto iter = std::sregex_iterator(line.begin(), line.end(), space_delimited_uint_regex);
        for (auto i=iter; i!=std::sregex_iterator(); ++i) ++out;
        assert(out);
        return out;
    }

    static size_t nreal_given(const std::string &filename, const std::regex &data_regex){
        FileIterator fi(filename, data_regex);
        std::string line{fi.next()};
        std::smatch match;
        std::regex_search(line, match, complex_regex);
        if constexpr(m_nreal==1){
            // we mustn't try to store a complex-valued tensor in a real container
            assert(!match.size());
        }
        if (match.size()) return 2;
        else return 1;
    }

    static size_t indsfirst(const std::string &filename){
        FileIterator fi(filename, float_inds_regex);
        std::string line{fi.next()};
        std::smatch match;
        std::regex_search(line, match, inds_first_regex);
        return match.size()>0;
    }

    static bool extract_line(const std::string &line, const size_t &nind,
                             const size_t &nreal_given, const bool &indsfirst, defs::inds &inds, T &value){
        assert(inds.size()>=nind);
        std::istringstream iss(line);
        std::vector<std::string> split(std::istream_iterator<std::string>{iss},
                                       std::istream_iterator<std::string>());
        for (auto i{0ul}; i<nind; ++i){
            inds[i] = atoi((split.begin()+(indsfirst?0:1)+i)->data())-1;
        }

        std::smatch match;
        std::regex_search(line, match, nreal_given==1?float_regex:complex_regex);
        if (!match.size()) return false;
        iss = std::istringstream(match.str());
        iss >> value;
        return true;
    }

    static defs::inds shape(const std::string &filename, const size_t &nind,
                            const size_t &nreal_given, const bool &indsfirst, const std::regex &data_regex){
        FileIterator fi(filename, data_regex);
        std::string line;
        defs::inds shape(nind, 0ul);
        defs::inds inds(nind, 0ul);
        T value;
        while (fi.next(line)){
            extract_line(line, nind, nreal_given, indsfirst, inds, value);
            for (auto i{0ul}; i<nind; ++i) {
                if (inds[i]<(size_t)(-1) && inds[i]+1>shape[i]) shape[i]=inds[i]+1;
            }
        }
        return shape;
    }

};

#endif //M7_TENSORFILEITERATOR_H
