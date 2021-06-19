//
// Created by Robert John Anderson on 2020-01-04.
//

#ifndef M7_LOGGING_H
#define M7_LOGGING_H

#include <iostream>
#include <execinfo.h>
#include "src/core/parallel/MPIWrapper.h"
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

    static void initialize(){
        if (!mpi::initialized())
            throw std::runtime_error("Logging requires that the MPIWrapper is initialized");
        if (mpi::i_am_root()) {
            g_reduced_stdout_logger = spdlog::stdout_color_st("stdout");
            g_reduced_file_logger = spdlog::basic_logger_st("fileout", "M7.log", true);
        }
#ifdef ENABLE_LOCAL_LOGGING
        g_local_file_logger = spdlog::basic_logger_st("fileout_", "M7.log."+std::to_string(mpi::irank()), true);
#endif
        spdlog::set_pattern("[%E] %^[%l]%$ %v");

#ifndef NDEBUG
        spdlog::set_level(spdlog::level::debug);
#else
        spdlog::set_level(spdlog::level::info);
#endif
#ifdef ENABLE_LOCAL_LOGGING
        g_local_file_logger->flush_on(spdlog::level::debug);
#endif
    }

    static void flush(){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->flush();
        g_reduced_file_logger->flush();
    }

    static void flush_(){
#ifdef ENABLE_LOCAL_LOGGING
        g_local_file_logger->flush();
#endif
    }

    static void flush_all(){
        flush();
        flush_();
    }

    static void finalize(){
        // spdlog::shutdown() doesn't seem to flush all streams in synchronous mode
        flush_all();
        spdlog::shutdown();
    }


    template<typename ...Args>
    static std::string format(const std::string& fmt_string, Args... args){
        fmt::basic_memory_buffer<char, 250> buf;
        fmt::format_to(buf, fmt_string, std::forward<Args>(args)...);
        std::string tmp;
        tmp.assign(buf.data(), buf.size());
        return tmp;
    }

    static std::vector<std::string> get_backtrace(size_t depth){
        std::vector<void*> entries(depth);
        size_t size;
        size = backtrace(entries.data(), depth);
        auto symbols = backtrace_symbols(entries.data(), size);
        std::vector<std::string> tmp;
        for (size_t i=0ul; i<size; ++i) tmp.emplace_back(symbols[i]);
        return tmp;
    }

    template<typename ...Args>
    static void info(const std::string& fmt_string, Args... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->info(fmt_string, args...);
        g_reduced_file_logger->info(fmt_string, args...);
    }

    template<typename ...Args>
    static void info_(const std::string& fmt_string, Args... args){
#ifdef ENABLE_LOCAL_LOGGING
        g_local_file_logger->info(fmt_string, args...);
#endif
    }

    template<typename ...Args>
    static void warn(const std::string& fmt_string, Args... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->warn(fmt_string, args...);
        g_reduced_file_logger->warn(fmt_string, args...);
    }

    template<typename ...Args>
    static void warn_(const std::string& fmt_string, Args... args){
#ifdef ENABLE_LOCAL_LOGGING
        g_local_file_logger->warn(fmt_string, args...);
#endif
    }

    template<typename ...Args>
    static void error(const std::string& fmt_string, Args... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->error(fmt_string, args...);
        g_reduced_file_logger->error(fmt_string, args...);
    }

    template<typename ...Args>
    static void error_(const std::string& fmt_string, Args... args){
#ifdef ENABLE_LOCAL_LOGGING
        g_local_file_logger->error(fmt_string, args...);
#endif
    }

    static void error_backtrace_(size_t depth=10){
        auto tmp = get_backtrace(depth);
        std::string str;
        error_("backtrace:");
        for (const auto& line: tmp) error_(line);
    }

    template<typename ...Args>
    static void critical(const std::string& fmt_string, Args... args){
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->critical(fmt_string, args...);
        g_reduced_file_logger->critical(fmt_string, args...);
    }

    template<typename ...Args>
    static void critical_(const std::string& fmt_string, Args... args){
#ifdef ENABLE_LOCAL_LOGGING
        g_local_file_logger->critical(fmt_string, args...);
#endif
    }

    template<typename ...Args>
    static void debug(const std::string& fmt_string, Args... args){
#ifndef NDEBUG
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->debug(fmt_string, args...);
        g_reduced_file_logger->debug(fmt_string, args...);
#endif
    }

    template<typename ...Args>
    static void debug_(const std::string& fmt_string, Args... args){
#ifndef NDEBUG
#ifdef ENABLE_LOCAL_LOGGING
        g_local_file_logger->debug(fmt_string, args...);
#endif
#endif
    }
};

#endif