//
// Created by Robert J. Anderson on 25/06/2021.
//

#include "YamlWrapper.h"
#include "M7_lib/util/String.h"

yaml::Path::Path(std::list<std::string> name_list) : m_name_list(std::move(name_list)) {}

yaml::Path::Path(std::vector<std::string> name_list) : m_name_list(name_list.cbegin(), name_list.cend()) {}

yaml::Path::Path(std::string name) : Path(string::split(name, '.')) {}

yaml::Path::Path(const yaml::Path &other) : m_name_list(other.m_name_list) {}

std::string yaml::Path::to_string() const {
    if (m_name_list.empty()) return "";
    auto it = m_name_list.cbegin();
    std::string out = *it;
    ++it;
    for (; it != m_name_list.cend(); ++it) out += "." + *it;
    return out;
}

yaml::Path yaml::Path::operator+(const std::string &name) const {
    auto tmp = *this;
    tmp.m_name_list.push_back(name);
    return tmp;
}

yaml::Path yaml::Path::up() const {
    auto tmp = *this;
    tmp.m_name_list.pop_back();
    return tmp;
}

size_t yaml::Path::depth() const {
    return m_name_list.size();
}

yaml::File::File(const std::string &fname) : m_fname(fname) {
    std::string contents;
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

YAML::Node yaml::File::get(std::list<std::string> path) const {
    return get(Path{path});
}

YAML::Node yaml::File::get(std::string path) const {
    return get(Path{path});
}

bool yaml::File::exists(yaml::Path path) const {
    return get(path).IsDefined();
}

bool yaml::File::exists(std::list<std::string> path) const {
    return exists(Path{path});
}

bool yaml::File::exists(std::string path) const {
    return exists(Path{path});
}
