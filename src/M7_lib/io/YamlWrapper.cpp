//
// Created by Robert J. Anderson on 25/06/2021.
//

#include "YamlWrapper.h"
#include "M7_lib/util/String.h"

yaml::Path::Path(std::list<str_t> name_list) : m_name_list(std::move(name_list)) {}

yaml::Path::Path(strv_t name_list) : m_name_list(name_list.cbegin(), name_list.cend()) {}

yaml::Path::Path(str_t name) : Path(string::split(name, '.')) {}

yaml::Path::Path(const yaml::Path &other) : m_name_list(other.m_name_list) {}

str_t yaml::Path::to_string() const {
    if (m_name_list.empty()) return "";
    auto it = m_name_list.cbegin();
    str_t out = *it;
    ++it;
    for (; it != m_name_list.cend(); ++it) out += "." + *it;
    return out;
}

yaml::Path yaml::Path::operator+(const str_t &name) const {
    auto tmp = *this;
    tmp.m_name_list.push_back(name);
    return tmp;
}

yaml::Path yaml::Path::up() const {
    auto tmp = *this;
    tmp.m_name_list.pop_back();
    return tmp;
}

uint_t yaml::Path::depth() const {
    return m_name_list.size();
}

yaml::File::File(const str_t &fname) : m_fname(fname) {
    str_t contents;
    /*
     * read in YAML file contents on the root node then then bcast to others
     */
    if (mpi::i_am_root()) contents = FileReader::to_string(fname);
    mpi::bcast(contents);
    try {
        m_root = YAML::Load(contents);
    }
    catch (const YAML::ParserException& ex) {
        ABORT(log::format("YAML syntax error in file {}, line {}, column {}",
                          m_fname, ex.mark.line, ex.mark.pos));
    }
}

YAML::Node yaml::File::get(yaml::Path path) const {
    auto node = m_root;
    for (const auto &it : path.m_name_list) node.reset(node[it]);
    return node;
}

YAML::Node yaml::File::get(std::list<str_t> path) const {
    return get(Path{path});
}

YAML::Node yaml::File::get(str_t path) const {
    return get(Path{path});
}

bool yaml::File::exists(yaml::Path path) const {
    return get(path).IsDefined();
}

bool yaml::File::exists(std::list<str_t> path) const {
    return exists(Path{path});
}

bool yaml::File::exists(str_t path) const {
    return exists(Path{path});
}
