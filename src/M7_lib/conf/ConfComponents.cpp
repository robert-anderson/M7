//
// Created by Robert J. Anderson on 25/06/2021.
//

#include "ConfComponents.h"
#include "M7_lib/util/String.h"

#if 0
yaml::Node::Node(yaml::Node *parent, str_t name, str_t description, bool impl_enable) :
        m_parent(parent), m_yaml_path(parent ? yaml::Path(parent->m_yaml_path + name) : yaml::Path(name)),
        m_description(description), m_indent(2 * m_yaml_path.depth(), ' '), m_impl_enable(impl_enable) {
}

yaml::Node::Node(str_t description, bool impl_enable) :
    Node(nullptr, "", description, impl_enable){}

str_t yaml::Node::help_string() const {
    return "";
}

const yaml::Document *yaml::Node::get_file() const {
    if (!m_parent) return nullptr;
    return m_parent->get_file();
}

bool yaml::Node::exists_in_file() const {
    auto file = get_file();
    return file->exists(m_yaml_path);
}

str_t yaml::Node::invalid_file_key() const {
    return "";
}

const str_t &yaml::Node::name() const {
    return m_yaml_path.m_name_list.back();
}

bool yaml::Node::enabled_internal() const {
    return true;
}

bool yaml::Node::enabled() const {
    /*
     * only need to check for existence in file for the current node. if this node is in the file, then all its
     * ancestors must also be
     */
    if (!m_impl_enable && !exists_in_file()) return false;
    for (auto node = m_parent; node!= nullptr; node=node->m_parent){
        if (!node->enabled_internal()) return false;
    }
    return true;
}

std::set<str_t> yaml::Group::make_file_keys() const {
    std::set<str_t> file_keys;
    auto yf = get_file();
    if (!yf) return file_keys;
    auto yaml_node = yf->get(m_yaml_path);
    for (auto it = yaml_node.begin(); it != yaml_node.end(); ++it)
        file_keys.insert(YAML::Dump(it->first));
    return file_keys;
}

std::set<str_t> yaml::Group::make_child_keys() const {
    std::set<str_t> child_keys;
    for (auto child: m_children) {
        auto child_key = child->name();
        REQUIRE_FALSE(child_keys.count(child_key), "Shouldn't have two child Nodes with the same name");
        child_keys.insert(child_key);
    }
    return child_keys;
}

yaml::Group::Group(yaml::Group *parent, str_t name, str_t description, bool impl_enable) :
        Node(parent, name, description, impl_enable) {
    if (parent) parent->add_child(this);
}

yaml::Group::Group(str_t description) : Node(description, false) {}

void yaml::Group::add_child(yaml::Node *child) {
    m_children.push_back(child);
}

str_t yaml::Group::invalid_file_key() const {
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

yaml::Section::Section(yaml::Group *parent, str_t name, str_t description, bool impl_enabled) :
        Group(parent, name, description, impl_enabled) {}

str_t yaml::Section::help_string() const {
    str_t str;
    const str_t enable = m_impl_enable ? "implicit" : "explicit";
    str.append(log::format("{}Section:      {}\n", m_indent, log::bold_format(m_yaml_path.to_string())));
    str.append(log::format("{}Enable:       {}\n", m_indent, enable));
    str.append(log::format("{}Description:  {}\n\n", m_indent, m_description));
    for (auto child: m_children) str.append(child->help_string() + "\n");
    return str;
}

yaml::Document::Document(const yaml::Document *file, str_t name, str_t description) :
        Group(description), m_name(name), m_file(file) {}

const yaml::Document *yaml::Document::get_file() const {
    return m_file;
}

str_t yaml::Document::help_string() const {
    REQUIRE_FALSE(m_file, "Help string should only be generated when it has not been filled by a YAML file");
    str_t str = string::boxed(log::format("Configuration document \"{}\"", m_name));
    str.append(log::format("Description: {}\n\n", m_description));
    for (auto child: m_children) str.append(child->help_string());
    return str;
}

void yaml::Document::verify() {
    Group::verify();
    if (m_file) {
        auto invalid = invalid_file_key();
        REQUIRE_TRUE(invalid.empty(),
                     log::format("YAML file \"{}\" contains invalid key \"{}\"", m_file->m_fname, invalid));
    }
}

yaml::ParamBase::ParamBase(yaml::Group *parent, str_t name, str_t description,
                             str_t v_default_str, str_t dim_type_str) :
        Node(parent, name, description, false), m_v_default_str(v_default_str),
        m_dim_type_str(dim_type_str) {
    REQUIRE_TRUE(m_parent, "Non-root yaml::Nodes must have a parent");
    parent->add_child(this);
}

str_t yaml::ParamBase::help_string() const {
    str_t str;
    str.append(log::format("{}Parameter:        {}\n", m_indent, log::bold_format(name())));
    str.append(log::format("{}Type:             {}\n", m_indent, m_dim_type_str));
    str.append(log::format("{}Default value:    {}\n", m_indent, m_v_default_str));
    str.append(log::format("{}Description:      {}\n", m_indent, m_description));
    return str;
}

#endif