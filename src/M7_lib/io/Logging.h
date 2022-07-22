//
// Created by Robert John Anderson on 2020-01-04.
//

#ifndef M7_LOGGING_H
#define M7_LOGGING_H

#include <iostream>
#include <execinfo.h>
#include <cxxabi.h>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/fmt/ostr.h>

#include <M7_lib/parallel/MPIWrapper.h>
#include <M7_lib/util/Datatype.h>


/*
 * Two different kinds of log:
 *   1. Reduced: only output by the root MPI rank to stdout, and a non-color copy to M7.log
 *   2. Local: output by all ranks (including the root) to M7.log.<irank>
 */

extern std::shared_ptr<spdlog::logger> g_reduced_stdout_logger;
extern std::shared_ptr<spdlog::logger> g_reduced_file_logger;
extern std::shared_ptr<spdlog::logger> g_local_file_logger;

namespace logging {

    void initialize();

    void flush();

    void flush_();

    void flush_all();

    void finalize();

    template<typename ...Args>
    static str_t format(const str_t& fmt_string, Args&&... args){
        fmt::basic_memory_buffer<char, 250> buf;
        fmt::format_to(buf, fmt_string, std::forward<Args>(args)...);
        str_t tmp;
        tmp.assign(buf.data(), buf.size());
        return tmp;
    }

    template<typename ...Args>
    static str_t bold_format(const str_t& fmt_string, Args&&... args){
        return format("\033[1m{}\033[0m", format(fmt_string, args...));
    }

    str_t get_demangled_symbol(const str_t& symbol);

    str_t get_demangled_prototype(const char* line);

    strv_t get_backtrace(uint_t depth);

    template<typename ...Args>
    static void info(const str_t& fmt_string, Args&&... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->info(fmt_string, args...);
        g_reduced_file_logger->info(fmt_string, args...);
    }

    template<typename ...Args>
    static void info_(const str_t& fmt_string, Args&&... args){
        if (mpi::nrank() == 1) info(fmt_string, args...);
        if (c_enable_local_logging) {
            g_local_file_logger->info(fmt_string, args...);
        }
    }

    void info_lines(const strv_t& lines);

    void info_lines_(const strv_t& lines);

    /**
     * make a pretty table
     * @param rows
     *  vectors of string vectors, each string being a cell in the table. Rows do not have to be the same length
     * @param header
     *  if true, the initial row of the table is treated as a header and a special underline is applied
     * @param padding
     *  minimum whitespace between the text of a cell and the vertical divider character '|'
     * @return
     *  a vector of strings which when printed will display a table in which the rows are vertically aligned
     */
    strv_t make_table(const v_t<strv_t>& rows, bool header=false, uint_t padding=2);

    strv_t make_table(const str_t& title, const v_t<strv_t>& rows, bool header=false, uint_t padding=2);

    void info_table(const v_t<strv_t>& rows, bool header=false, uint_t padding=2);

    void info_table(const str_t& title, const v_t<strv_t>& rows, bool header=false, uint_t padding=2);

    void info_table_(const v_t<strv_t>& rows, bool header=false, uint_t padding=2);

    void info_table_(const str_t& title, const v_t<strv_t>& rows, bool header=false, uint_t padding=2);

    template<typename ...Args>
    static void warn(const str_t& fmt_string, Args&&... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->warn(fmt_string, args...);
        g_reduced_file_logger->warn(fmt_string, args...);
    }

    template<typename ...Args>
    static void warn_(const str_t& fmt_string, Args&&... args){
        if (mpi::nrank() == 1) warn(fmt_string, args...);
        if (c_enable_local_logging) {
            g_local_file_logger->warn(fmt_string, args...);
        }
    }

    template<typename ...Args>
    static void error(const str_t& fmt_string, Args&&... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->error(fmt_string, args...);
        g_reduced_file_logger->error(fmt_string, args...);
    }

    template<typename ...Args>
    static void error_(const str_t& fmt_string, Args&&... args){
        if (mpi::nrank() == 1) error(fmt_string, args...);
        if (c_enable_local_logging) {
            g_local_file_logger->error(fmt_string, args...);
        }
    }

    void error_backtrace_(uint_t depth=20);

    template<typename ...Args>
    static void critical(const str_t& fmt_string, Args&&... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->critical(fmt_string, args...);
        g_reduced_file_logger->critical(fmt_string, args...);
    }

    template<typename ...Args>
    static void critical_(const str_t& fmt_string, Args&&... args){
        if (mpi::nrank() == 1) critical(fmt_string, args...);
        if (c_enable_local_logging) {
            g_local_file_logger->critical(fmt_string, args...);
        }
    }

    template<typename ...Args>
    static void debug(const str_t& fmt_string, Args&&... args){
        if (c_enable_debug) {
            if (!mpi::i_am_root()) return;
            g_reduced_stdout_logger->debug(fmt_string, args...);
            g_reduced_file_logger->debug(fmt_string, args...);
        }
    }

    template<typename ...Args>
    static void debug_(const str_t& fmt_string, Args&&... args){
        if (mpi::nrank() == 1) debug(fmt_string, args...);
        if (c_enable_debug && c_enable_local_logging) {
            g_local_file_logger->debug(fmt_string, args...);
        }
    }

    str_t quoted(const str_t& str);

    template<uint_t mbf_ind=c_mbf_type_ind>
    static str_t mbf_name(){
        switch (mbf_ind) {
            case 0: return "fermion (determinant)";
            case 1: return "fermion-boson (determinant-permanent product)";
            case 2: return "boson (permanent)";
            case 3: return "fermion spin-adapted (CSF)";
            default: return "";
        }
    }

    /**
     * print the M7 title image
     */
    void title();

    /**
     * print the compile-time definitions to standard output
     */
    void defs();
}

#endif
