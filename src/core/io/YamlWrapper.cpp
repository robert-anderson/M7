//
// Created by rja on 25/06/2021.
//

#include "YamlWrapper.h"

yaml::Path::Path(std::list<std::string> path) : m_path(std::move(path)) {}

yaml::Path::Path(std::vector<std::string> path) : m_path(path.cbegin(), path.cend()) {}

yaml::Path::Path(std::string path) : Path(string_utils::split(path, '.')) {}

yaml::Path::Path(const yaml::Path &other) : m_path(other.m_path) {}

std::string yaml::Path::to_string() const {
    auto it = m_path.cbegin();
    std::string out = *it;
    ++it;
    for (; it != m_path.cend(); ++it) out += "." + *it;
    return out;
}

yaml::Path yaml::Path::operator+(const std::string &name) const {
    auto tmp = *this;
    tmp.m_path.push_back(name);
    return tmp;
}

yaml::Path yaml::Path::up() const {
    auto tmp = *this;
    tmp.m_path.pop_back();
    return tmp;
}

size_t yaml::Path::depth() const {
    return m_path.size();
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
        ABORT(log::format("YAML syntax error in file {}", m_fname));
    }
}

YAML::Node yaml::File::get(yaml::Path path) const {
    auto node = m_root;
    for (const auto &it : path.m_path) node.reset(node[it]);
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
