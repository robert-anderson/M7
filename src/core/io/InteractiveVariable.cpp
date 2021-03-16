//
// Created by rja on 16/03/2021.
//

#include "InteractiveVariable.h"


InteractiveVariableFile::InteractiveVariableFile(std::string name) : m_name(name), m_fname(name + ".var") {
    log::info(R"(Listening for changes to interactive variable "{}" in file "{}")", name, m_fname);
}

bool InteractiveVariableFile::consume_(std::vector<std::string> &lines) {
    if (mpi::i_am_root()) {
        std::ifstream f(m_fname);
        if (!f.is_open()) return false;
        do {
            lines.emplace_back();
            std::getline(f, lines.back());
        } while (!lines.back().empty());
        f.close();
        std::remove(m_fname.c_str());
        lines.erase(lines.end());
    }
    return true;
}

void InteractiveVariableFile::warn_invalid_input() const {
    log::warn("Invalid value given in file \"{}\", variable left unchanged", m_fname);
}

void InteractiveVariableFile::info_success(const std::string &str) const {
    log::info("Interactive variable \"{}\" assigned value: {}", m_name, str);
}
