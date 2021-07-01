//
// Created by rja on 26/01/2021.
//

#include "Logging.h"

std::shared_ptr<spdlog::logger> g_reduced_stdout_logger = nullptr;
std::shared_ptr<spdlog::logger> g_reduced_file_logger = nullptr;
std::shared_ptr<spdlog::logger> g_local_file_logger = nullptr;

void log::initialize() {
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

void log::flush() {
    if (!mpi::i_am_root()) return;
    g_reduced_stdout_logger->flush();
    g_reduced_file_logger->flush();
}

void log::flush_() {
#ifdef ENABLE_LOCAL_LOGGING
    g_local_file_logger->flush();
#endif
}

void log::flush_all() {
    flush();
    flush_();
}

void log::finalize() {
    // spdlog::shutdown() doesn't seem to flush all streams in synchronous mode
    flush_all();
    spdlog::shutdown();
}

std::string log::get_demangled_symbol(const std::string &symbol) {
    int status;
    auto demangled = abi::__cxa_demangle(symbol.data(), nullptr, nullptr, &status);
    if (status!=0) {
        free(demangled);
        return "";
    }
    std::string tmp(demangled);
    free(demangled);
    return tmp;
}

std::string log::get_demangled_prototype(const char *line) {
    // parse till first "(", then till "+"
    const char* begin_ptr = nullptr;
    const char* end_ptr = nullptr;
    for (const char* ptr=line; *ptr!=0; ++ptr){
        if (!begin_ptr && *ptr=='(') begin_ptr = ptr+1;
        else if (begin_ptr && *ptr=='+') {
            end_ptr = ptr;
            break;
        }
    }
    size_t length = std::distance(begin_ptr, end_ptr);
    std::string symbol(begin_ptr, length);
    return get_demangled_symbol(symbol);
}

std::vector<std::string> log::get_backtrace(size_t depth) {
    std::vector<void*> entries(depth);
    size_t size;
    size = backtrace(entries.data(), depth);
    auto symbols = backtrace_symbols(entries.data(), size);
    std::vector<std::string> tmp;
    for (size_t i=0ul; i<size; ++i) {
        tmp.emplace_back(get_demangled_prototype(symbols[i]));
    }
    free(symbols);
    return tmp;
}

void log::error_backtrace_(size_t depth) {
    auto tmp = get_backtrace(depth);
    std::string str;
    error_("Printing backtrace (maximum call depth {}):", depth);
    for (const auto& line: tmp) if (!line.empty()) error_(line);
}
