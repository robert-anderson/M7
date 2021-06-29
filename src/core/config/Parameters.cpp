//
// Created by rja on 25/06/2021.
//

#include "Parameters.h"

config::Node::Node(config::Node *parent, std::string name, std::string description) :
        m_parent(parent), m_path(parent->m_path + name),
        m_description(description), m_indent(2 * (m_path.depth() - 1), ' ') {
}

config::Node::Node(std::string description) :
        m_parent(nullptr), m_path(""), m_description(description), m_indent() {}

std::string config::Node::help_string() const {
    return "";
}

const yaml::File *config::Node::get_file() const {
    REQUIRE_TRUE(m_parent, "Non-root Nodes must have a parent");
    return m_parent->get_file();
}

std::string config::Node::invalid_file_key() const {
    return "";
}

const std::string &config::Node::name() const {
    return m_path.m_path.back();
}

std::set<std::string> config::Group::make_file_keys() const {
    std::set<std::string> file_keys;
    auto yf = get_file();
    if (!yf) return file_keys;
    auto yaml_node = yf->get(m_path);
    for (auto it = yaml_node.begin(); it != yaml_node.end(); ++it)
        file_keys.insert(YAML::Dump(it->first));
    return file_keys;
}

std::set<std::string> config::Group::make_child_keys() const {
    std::set<std::string> child_keys;
    for (auto child: m_children) {
        auto child_key = child->name();
        REQUIRE_FALSE(child_keys.count(child_key), "Shouldn't have two child Nodes with the same name");
        child_keys.insert(child_key);
    }
    return child_keys;
}

config::Group::Group(config::Group *parent, std::string name, std::string description) :
        Node(parent, name, description) {
    REQUIRE_TRUE(m_parent, "Non-root config::Nodes must have a parent");
    parent->add_child(this);
}

config::Group::Group(std::string description) : Node(description) {}

void config::Group::add_child(config::Node *child) {
    m_children.push_back(child);
}

std::string config::Group::invalid_file_key() const {
    if (get_file()) {
        auto file_keys = make_file_keys();
        auto child_keys = make_child_keys();
        for (const auto &it : file_keys)
            if (!child_keys.count(it)) return it;
    }
    for (auto child : m_children) {
        auto invalid = child->invalid_file_key();
        if (!invalid.empty()) return invalid;
    }
    return "";
}

bool config::Section::make_exists() const {
    auto file = get_file();
    if (file) return file->exists(m_path);
    return false;
}

config::Section::Section(config::Group *parent, std::string name, std::string description) :
        Group(parent, name, description), m_is_specified(make_exists()) {}

config::Section::operator bool() const {
    return m_is_specified;
}

std::string config::Section::help_string() const {
    std::string str;
    str.append(log::format("{}Section:       {}\n", m_indent, log::bold_format(m_path.to_string())));
    str.append(log::format("{}Description:   {}\n\n", m_indent, m_description));
    for (auto child: m_children) str.append(child->help_string() + "\n");
    return str;
}

config::Document::Document(const yaml::File *file, std::string name, std::string description) :
        Group(description), m_name(name), m_file(file) {}

const yaml::File *config::Document::get_file() const {
    return m_file;
}

std::string config::Document::help_string() const {
    REQUIRE_FALSE(m_file, "Help string should only be generated when it has not been filled by a YAML file");
    std::string str = string_utils::boxed(log::format("Configuration document \"{}\"", m_name));
    str.append(log::format("Description: {}\n\n", m_description));
    for (auto child: m_children) str.append(child->help_string());
    return str;
}

void config::Document::verify() {
    Node::verify();
    if (m_file) {
        auto invalid = invalid_file_key();
        REQUIRE_TRUE(invalid.empty(),
                     log::format("YAML file \"{}\" contains invalid key \"{}\"", m_file->m_fname, invalid));
    }
}

config::ParamBase::ParamBase(config::Group *parent, std::string name, std::string description,
                             std::string v_default_str, std::string dim_type_str) :
        Node(parent, name, description), m_v_default_str(v_default_str),
        m_dim_type_str(dim_type_str) {
    REQUIRE_TRUE(m_parent, "Non-root config::Nodes must have a parent");
    parent->add_child(this);
}

std::string config::ParamBase::help_string() const {
    std::string str;
    str.append(log::format("{}Parameter:       {}\n", m_indent, log::bold_format(name())));
    str.append(log::format("{}Type:            {}\n", m_indent, m_dim_type_str));
    str.append(log::format("{}Default value:   {}\n", m_indent, m_v_default_str));
    str.append(log::format("{}Description:     {}\n", m_indent, m_description));
    return str;
}
