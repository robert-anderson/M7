//
// Created by Robert J. Anderson on 26/01/2021.
//

#include "Logging.h"
#include "M7_lib/version/version.h"

std::shared_ptr<spdlog::logger> g_reduced_stdout_logger = nullptr;
std::shared_ptr<spdlog::logger> g_reduced_file_logger = nullptr;
std::shared_ptr<spdlog::logger> g_local_file_logger = nullptr;

void logging::initialize() {
    if (!mpi::initialized())
        throw std::runtime_error("Logging requires that the MPIWrapper is initialized");
    if (mpi::i_am_root()) {
        g_reduced_stdout_logger = spdlog::stdout_color_st("stdout");
        g_reduced_file_logger = spdlog::basic_logger_st("fileout", "M7.log", true);
    }
    if (c_enable_local_logging)
        g_local_file_logger = spdlog::basic_logger_st("fileout_", "M7.log."+std::to_string(mpi::irank()), true);
    spdlog::set_pattern("[%E] %^[%l]%$ %v");

    spdlog::set_level(c_enable_debug ? spdlog::level::debug : spdlog::level::info);
    if (c_enable_local_logging) g_local_file_logger->flush_on(spdlog::level::debug);
}

void logging::flush() {
    if (!mpi::i_am_root()) return;
    g_reduced_stdout_logger->flush();
    g_reduced_file_logger->flush();
}

void logging::flush_() {
    if (c_enable_local_logging) g_local_file_logger->flush();
}

void logging::flush_all() {
    flush();
    flush_();
}

void logging::finalize() {
    // spdlog::shutdown() doesn't seem to flush all streams in synchronous mode
    flush_all();
    spdlog::shutdown();
}

str_t logging::get_demangled_symbol(const str_t &symbol) {
    int status;
    auto demangled = abi::__cxa_demangle(symbol.data(), nullptr, nullptr, &status);
    if (status!=0) {
        free(demangled);
        return "";
    }
    str_t tmp(demangled);
    free(demangled);
    return tmp;
}

str_t logging::get_demangled_prototype(const char *line) {
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
    uint_t length = std::distance(begin_ptr, end_ptr);
    str_t symbol(begin_ptr, length);
    return get_demangled_symbol(symbol);
}

strv_t logging::get_backtrace(uint_t depth) {
    v_t<void*> entries(depth);
    uint_t size;
    size = backtrace(entries.data(), depth);
    auto symbols = backtrace_symbols(entries.data(), size);
    strv_t tmp;
    for (uint_t i=0ul; i<size; ++i) {
        tmp.emplace_back(get_demangled_prototype(symbols[i]));
    }
    free(symbols);
    return tmp;
}

void logging::error_backtrace_(uint_t depth) {
    auto tmp = get_backtrace(depth);
    str_t str;
    error_("Printing backtrace (maximum call depth {}):", depth);
    for (const auto& line: tmp) if (!line.empty()) error_(line);
}

strv_t logging::make_table(const v_t<strv_t> &rows, bool header, bool hlines, uint_t padding) {
    if (rows.empty()) return {};
    auto fn = [](const strv_t& row1, const strv_t& row2){
        return row1.size() < row2.size();
    };

    auto max_it = std::max_element(rows.cbegin(), rows.cend(), fn);
    const uint_t ncol = max_it->size();
    uintv_t max_sizes(ncol, 0ul);
    for (auto& row: rows) {
        for (auto it=row.cbegin(); it!=row.cend(); ++it) {
            auto icol = std::distance(row.cbegin(), it);
            max_sizes[icol] = std::max(max_sizes[icol], it->size());
        }
    }
    str_t hline{'+'};
    str_t header_hline{'+'};
    for (auto& max_size: max_sizes) {
        hline.append(max_size+2*padding,'-');
        hline.append(1, '+');
        if (header) {
            header_hline.append(max_size + 2 * padding, '=');
            header_hline.append(1, '+');
        }
    }
    strv_t row_strs;
    row_strs.reserve(2*rows.size()+1);
    row_strs.push_back(hline);
    for (auto& row: rows) {
        str_t row_str = "|";
        for (uint_t icol=0ul; icol<ncol; ++icol){
            row_str.append(padding, ' ');
            row_str.append(icol<row.size() ? row[icol] : "");
            row_str.append(padding + max_sizes[icol]-(icol<row.size() ? row[icol].size() : 0ul), ' ');
            row_str.append(1, '|');
        }
        row_strs.push_back(row_str);
        if (row_strs.size()==2 && header) row_strs.push_back(header_hline);
        else if (hlines) row_strs.push_back(hline);
    }
    // even if there are no hlines between rows, place one at the end
    if (!hlines) row_strs.push_back(hline);
    return row_strs;
}

strv_t logging::make_table(const str_t &title, const v_t<strv_t> &rows, bool header, bool hlines, uint_t padding) {
    if (rows.empty()) return {};
    strv_t table;
    auto body = make_table(rows, header, hlines, padding);
    table.push_back(body.front());
    str_t title_str = "|";
    const uint_t titlebar_size = body.front().size()-2;
    const uint_t nspace_left = titlebar_size < title.size() ? 0ul : (titlebar_size - title.size()) / 2;
    const uint_t nspace_right = titlebar_size < title.size() ? 0ul : titlebar_size - (title.size() + nspace_left);
    title_str.append(nspace_left, ' ');
    title_str.insert(title_str.end(), title.cbegin(), nspace_right ? title.cend() : title.cbegin() + titlebar_size);
    title_str.append(nspace_right, ' ');
    title_str.append(1, '|');
    table.push_back(title_str);
    table.insert(table.end(), body.cbegin(), body.cend());
    return table;
}

void logging::info_lines(const strv_t &lines) {
    for (const auto& line: lines) info(line);
}

void logging::info_lines_(const strv_t &lines) {
    for (const auto& line: lines) info_(line);
}

void logging::info_table(const v_t<strv_t> &rows, bool header, bool hlines, uint_t padding) {
    info_lines(make_table(rows, header, hlines, padding));
}

void logging::info_table(const str_t &title, const v_t<strv_t> &rows, bool header, bool hlines, uint_t padding) {
    info_lines(make_table(title, rows, header, hlines, padding));
}

void logging::info_table_(const v_t<strv_t> &rows, bool header, bool hlines, uint_t padding) {
    info_lines_(make_table(rows, header, hlines, padding));
}

void logging::info_table_(const str_t &title, const v_t<strv_t> &rows, bool header, bool hlines, uint_t padding) {
    info_lines_(make_table(title, rows, header, hlines, padding));
}

void logging::title() {
    if (!mpi::i_am_root()) return;
    std::cout << R"(   ____    ____   _________              )" << '\n';
    std::cout << R"(  |****\  /****| |********/  M any-body  )" << '\n';
    std::cout << R"(  |**|\*\/*/|**|      /**/   S tochastic )" << '\n';
    std::cout << R"(  |**| \**/ |**|     /**/    E xpectation)" << '\n';
    std::cout << R"(  |**|      |**|    /**/     V alue      )" << '\n';
    std::cout << R"(  |**|      |**|   /**/      E stimation )" << '\n';
    std::cout << R"(  |**|      |**|  /**/       N etworks   )" << '\n';
    std::cout << std::endl;
}

str_t logging::quoted(const str_t &str) {
    return "\""+str+"\"";
}

strv_t logging::make_defs_table() {
    return make_table("compile definitions",{
        {"git revision", get_version()},
        {"build mode", c_enable_debug ? "DEBUG" : "RELEASE"},
        {"compilation timestamp", get_compilation_timestamp()},
        {"many-body basis function", mbf_name<c_mbf_type_ind>()},
        {"walker arithmetic", c_enable_complex_wf ? "complex" : "real"},
        {"Hamiltonian arithmetic", c_enable_complex_ham ? "complex" : "real"},
        {"TCHINT interface", c_enable_tchint ? "enabled" : "disabled"},
        {"rank-resolved logging", c_enable_local_logging ? "enabled" : "disabled"},
    });
}
