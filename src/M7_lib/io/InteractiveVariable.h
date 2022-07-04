//
// Created by Robert J. Anderson on 16/03/2021.
//

#ifndef M7_INTERACTIVEVARIABLE_H
#define M7_INTERACTIVEVARIABLE_H

#include "FileReader.h"
#include "Logging.h"
#include "M7_lib/util/String.h"

/*
 * variables can be updated from a file:
 * if the interactive variable name is "shift", the file "shift.var" will be opened
 */

struct InteractiveVariableFile {
    str_t m_name, m_fname;

    InteractiveVariableFile(str_t name);

private:
    bool consume_(strv_t &lines);

    void warn_invalid_input() const;

    void info_success(const str_t &str) const;

    template<typename T>
    typename std::enable_if<std::is_integral<T>::value && std::is_signed<T>::value, bool>::type
    read_line(const str_t &line, T &v) {
        auto ptr = line.c_str();
        int64_t tmp = string::read_signed(ptr);
        v = static_cast<T>(tmp);
        return static_cast<const int64_t &>(v) == tmp;
    }

    template<typename T>
    typename std::enable_if<std::is_integral<T>::value && !std::is_signed<T>::value, bool>::type
    read_line(const str_t &line, T &v) {
        if (line[0] == '-') return false;
        auto ptr = line.c_str();
        uint_t tmp = string::read_unsigned(ptr);
        v = static_cast<T>(tmp);
        return static_cast<const uint_t &>(v) == tmp;
    }


    template<typename T>
    typename std::enable_if<std::is_floating_point<T>::value, bool>::type
    read_line(const str_t &line, T &v) {
        auto ptr = line.c_str();
        auto tmp = string::read_double(ptr);
        if (tmp==std::numeric_limits<double>::max()) return false;
        v = tmp;
        return true;
    }

    template<typename T>
    typename std::enable_if<std::is_floating_point<T>::value, bool>::type
    read_line(const str_t &line, std::complex<T> &v) {
        auto ptr = line.c_str();
        auto real = string::read_double(ptr);
        if (real==std::numeric_limits<double>::max()) return false;
        auto imag = *ptr==0 ? 0.0: string::read_double(ptr);
        if (imag==std::numeric_limits<double>::max()) return false;
        v = {real, imag};
        return true;
    }

    template<typename T>
    bool read_vector(v_t<T> &v, strv_t& lines) {
        uint_t nelement = lines.size();
        mpi::bcast(nelement);
        v_t<T> tmp;
        tmp.reserve(nelement);
        bool invalid = false;
        if (mpi::i_am_root()) {
            for (const auto &line: lines) {
                tmp.push_back({});
                invalid = !read_line(line, tmp.back());
                if (invalid) break;
            }
        }
        mpi::bcast(invalid, 0);
        if (invalid) return false;
        v = tmp;
        mpi::bcast(v, nelement, 0);
        return true;
    }

    template<typename T>
    bool read_check(T &v, strv_t& lines) {
        v_t<T> tmp{v};
        auto res = read_vector(tmp, lines);
        if (!res) return false;
        v = tmp[0];
        return true;
    }

    bool read_check(v_t<bool> &v, strv_t& lines) {
        v_t<char> tmp;
        auto res = read_vector(tmp, lines);
        if (!res) return false;
        v_t<bool> tmp2;
        tmp.reserve(tmp.size());
        for (const auto &c: tmp) {
            if (c == 0 || c == 1) tmp2.push_back(!!c);
            else return false;
        }
        if (tmp2.size()==v.size()) v = tmp2;
        else if (tmp2.size()==1) v.assign(v.size(), tmp[0]);
        else return false;
        return true;
    }

    bool read_check(bool &v, strv_t& lines) {
        v_t<bool> tmp{v};
        auto res = read_check(tmp, lines);
        if (!res) return false;
        v = tmp[0];
        return true;
    }

    template<typename T>
    bool read_check(v_t<T> &v, strv_t& lines) {
        v_t<T> tmp{v};
        auto res = read_vector(tmp, lines);
        if (!res) return false;
        if (tmp.size()==v.size()) v = tmp;
        else if (tmp.size()==1) v.assign(v.size(), tmp[0]);
        else return false;
        return true;
    }

public:
    /**
     * @tparam T
     * type of variable being updated
     * @param v
     * ref to variable being updated
     * @return
     * true if the variable was updated, i.e. the file existed and contained valid values
     * else false
     */
    template<typename T>
    bool read(T &v) {
        strv_t lines;
        bool exists = consume_(lines);
        mpi::bcast(exists, 0);
        if (!exists) return false;
        bool res = read_check(v, lines);
        if (!res) {
            warn_invalid_input();
            return false;
        }
        else {
            info_success(convert::to_string(v));
            return true;
        }
    }
};

template<typename T>
struct InteractiveVariable {
    InteractiveVariableFile m_ivf;
    T m_v;
    InteractiveVariable(str_t name): m_ivf(name){}
    InteractiveVariable(str_t name, const T& v): InteractiveVariable(name){
        m_v = v;
    }

    operator T&() {
        return m_v;
    }

    operator const T&() const {
        return m_v;
    }

    bool read(){
        return m_ivf.read(m_v);
    }
};

#endif //M7_INTERACTIVEVARIABLE_H
