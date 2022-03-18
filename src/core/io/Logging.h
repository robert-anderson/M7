//
// Created by Robert John Anderson on 2020-01-04.
//

#ifndef M7_LOGGING_H
#define M7_LOGGING_H

#include <iostream>
#include <execinfo.h>
#include <cxxabi.h>
#include "parallel/MPIWrapper.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/fmt/ostr.h"


/*
 * Two different kinds of log:
 *   1. Reduced: only output by the root MPI rank to stdout, and a non-color copy to M7.log
 *   2. Local: output by all ranks (including the root) to M7.log.<irank>
 */

extern std::shared_ptr<spdlog::logger> g_reduced_stdout_logger;
extern std::shared_ptr<spdlog::logger> g_reduced_file_logger;
extern std::shared_ptr<spdlog::logger> g_local_file_logger;

struct log {

    static void initialize();

    static void flush();

    static void flush_();

    static void flush_all();

    static void finalize();

    template<typename ...Args>
    static std::string format(const std::string& fmt_string, Args&&... args){
        fmt::basic_memory_buffer<char, 250> buf;
        fmt::format_to(buf, fmt_string, std::forward<Args>(args)...);
        std::string tmp;
        tmp.assign(buf.data(), buf.size());
        return tmp;
    }

    template<typename ...Args>
    static std::string bold_format(const std::string& fmt_string, Args&&... args){
        return format("\033[1m{}\033[0m", format(fmt_string, args...));
    }

    static std::string get_demangled_symbol(const std::string& symbol);

    static std::string get_demangled_prototype(const char* line);

    static std::vector<std::string> get_backtrace(size_t depth);

    template<typename ...Args>
    static void info(const std::string& fmt_string, Args&&... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->info(fmt_string, args...);
        g_reduced_file_logger->info(fmt_string, args...);
    }

    template<typename ...Args>
    static void info_(const std::string& fmt_string, Args&&... args){
#ifdef ENABLE_LOCAL_LOGGING
        if (mpi::nrank()==1) info(fmt_string, args...);
        g_local_file_logger->info(fmt_string, args...);
#endif
    }

    template<typename ...Args>
    static void warn(const std::string& fmt_string, Args&&... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->warn(fmt_string, args...);
        g_reduced_file_logger->warn(fmt_string, args...);
    }

    template<typename ...Args>
    static void warn_(const std::string& fmt_string, Args&&... args){
#ifdef ENABLE_LOCAL_LOGGING
        if (mpi::nrank()==1) warn(fmt_string, args...);
        g_local_file_logger->warn(fmt_string, args...);
#endif
    }

    template<typename ...Args>
    static void error(const std::string& fmt_string, Args&&... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->error(fmt_string, args...);
        g_reduced_file_logger->error(fmt_string, args...);
    }

    template<typename ...Args>
    static void error_(const std::string& fmt_string, Args&&... args){
#ifdef ENABLE_LOCAL_LOGGING
        if (mpi::nrank()==1) error(fmt_string, args...);
        g_local_file_logger->error(fmt_string, args...);
#endif
    }

    static void error_backtrace_(size_t depth=20);

    template<typename ...Args>
    static void critical(const std::string& fmt_string, Args&&... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->critical(fmt_string, args...);
        g_reduced_file_logger->critical(fmt_string, args...);
    }

    template<typename ...Args>
    static void critical_(const std::string& fmt_string, Args&&... args){
#ifdef ENABLE_LOCAL_LOGGING
        if (mpi::nrank()==1) critical(fmt_string, args...);
        g_local_file_logger->critical(fmt_string, args...);
#endif
    }

    template<typename ...Args>
    static void debug(const std::string& fmt_string, Args&&... args){
#ifndef NDEBUG
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->debug(fmt_string, args...);
        g_reduced_file_logger->debug(fmt_string, args...);
#endif
    }

    template<typename ...Args>
    static void debug_(const std::string& fmt_string, Args&&... args){
#ifndef NDEBUG
#ifdef ENABLE_LOCAL_LOGGING
        if (mpi::nrank()==1) debug(fmt_string, args...);
        g_local_file_logger->debug(fmt_string, args...);
#endif
#endif
    }
};

#endif