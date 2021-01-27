//
// Created by Robert John Anderson on 2020-01-04.
//

#ifndef M7_LOGGING_H
#define M7_LOGGING_H

#include <iostream>
#include "src/core/parallel/MPIWrapper.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"


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
        spdlog::flush_every(std::chrono::seconds(3));
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
#ifndef DNDEBUG
        if (!mpi::i_am_root()) return;
        g_reduced_stdout_logger->debug(fmt_string, args...);
        g_reduced_file_logger->debug(fmt_string, args...);
#endif
    }

    template<typename ...Args>
    static void debug_(const std::string& fmt_string, Args... args){
#ifndef DNDEBUG
#ifdef ENABLE_LOCAL_LOGGING
        g_local_file_logger->debug(fmt_string, args...);
#endif
#endif
    }
};

#endif