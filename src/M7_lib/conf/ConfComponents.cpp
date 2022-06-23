//
// Created by Robert J. Anderson on 25/06/2021.
//

#include "ConfComponents.h"
#include "M7_lib/util/String.h"

conf_components::Node::Node(conf_components::Node *parent, std::string name, std::string description) :
        m_parent(parent), m_yaml_path(parent ? yaml::Path(parent->m_yaml_path + name) : yaml::Path(name)),
        m_description(description), m_indent(2 * m_yaml_path.depth(), ' ') {
}

conf_components::Node::Node(std::string description) : Node(nullptr, "", description){}

std::string conf_components::Node::help_string() const {
    return "";
}

const yaml::File *conf_components::Node::get_file() const {
    if (!m_parent) return nullptr;
    return m_parent->get_file();
}

std::string conf_components::Node::invalid_file_key() const {
    return "";
}

const std::string &conf_components::Node::name() const {
    return m_yaml_path.m_name_list.back();
}

bool conf_components::Node::parents_enabled() const {
    for (auto node = m_parent; node!= nullptr; node=node->m_parent){
        if (!node->enabled()) return false;
    }
    return true;
}

std::set<std::string> conf_components::Group::make_file_keys() const {
    std::set<std::string> file_keys;
    auto yf = get_file();
    if (!yf) return file_keys;
    auto yaml_node = yf->get(m_yaml_path);
    for (auto it = yaml_node.begin(); it != yaml_node.end(); ++it)
        file_keys.insert(YAML::Dump(it->first));
    return file_keys;
}

std::set<std::string> conf_components::Group::make_child_keys() const {
    std::set<std::string> child_keys;
    for (auto child: m_children) {
        auto child_key = child->name();
        REQUIRE_FALSE(child_keys.count(child_key), "Shouldn't have two child Nodes with the same name");
        child_keys.insert(child_key);
    }
    return child_keys;
}

conf_components::Group::Group(conf_components::Group *parent, std::string name, std::string description) :
        Node(parent, name, description) {
    if (parent) parent->add_child(this);
}

conf_components::Group::Group(std::string description) : Node(description) {}

void conf_components::Group::add_child(conf_components::Node *child) {
    m_children.push_back(child);
}

std::string conf_components::Group::invalid_file_key() const {
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

conf_components::Section::Section(conf_components::Group *parent, std::string name, std::string description) :
        Group(parent, name, description) {}

std::string conf_components::Section::help_string() const {
    std::string str;
    str.append(log::format("{}Section:       {}\n", m_indent, log::bold_format(m_yaml_path.to_string())));
    str.append(log::format("{}Description:   {}\n\n", m_indent, m_description));
    for (auto child: m_children) str.append(child->help_string() + "\n");
    return str;
}

conf_components::Document::Document(const yaml::File *file, std::string name, std::string description) :
        Group(description), m_name(name), m_file(file) {}

const yaml::File *conf_components::Document::get_file() const {
    return m_file;
}

std::string conf_components::Document::help_string() const {
    REQUIRE_FALSE(m_file, "Help string should only be generated when it has not been filled by a YAML file");
    std::string str = string::boxed(log::format("Configuration document \"{}\"", m_name));
    str.append(log::format("Description: {}\n\n", m_description));
    for (auto child: m_children) str.append(child->help_string());
    return str;
}

void conf_components::Document::verify() {
    Group::verify();
    if (m_file) {
        auto invalid = invalid_file_key();
        REQUIRE_TRUE(invalid.empty(),
                     log::format("YAML file \"{}\" contains invalid key \"{}\"", m_file->m_fname, invalid));
    }
}

conf_components::ParamBase::ParamBase(conf_components::Group *parent, std::string name, std::string description,
                             std::string v_default_str, std::string dim_type_str) :
        Node(parent, name, description), m_v_default_str(v_default_str),
        m_dim_type_str(dim_type_str) {
    REQUIRE_TRUE(m_parent, "Non-root conf_components::Nodes must have a parent");
    parent->add_child(this);
}

std::string conf_components::ParamBase::help_string() const {
    std::string str;
    str.append(log::format("{}Parameter:       {}\n", m_indent, log::bold_format(name())));
    str.append(log::format("{}Type:            {}\n", m_indent, m_dim_type_str));
    str.append(log::format("{}Default value:   {}\n", m_indent, m_v_default_str));
    str.append(log::format("{}Description:     {}\n", m_indent, m_description));
    return str;
}
