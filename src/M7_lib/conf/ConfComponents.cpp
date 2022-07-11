//
// Created by rja on 10/07/22.
//

#include "ConfComponents.h"

conf_components::Path::Path(std::list<str_t> list) :
        m_list(std::move(list)),
        m_string(string::join(v_t<str_t>(m_list.cbegin(), m_list.cend()), ".")){}

conf_components::Path conf_components::Path::operator+(const str_t& name) const {
    auto tmp = m_list;
    tmp.push_back(name);
    return {tmp};
}

v_t<str_t> conf_components::Node::child_names_in_file() const {
    if (!m_yaml_node.IsDefined()) return {};
    if (m_yaml_node.IsSequence()) return {};
    v_t<str_t> out;
    out.reserve(m_yaml_node.size());
    for (auto it=m_yaml_node.begin(); it!=m_yaml_node.end(); ++it) {
        out.push_back(it->first.Scalar());
    }
    return out;
}

void conf_components::Node::validate_file_contents() const {
    if (!m_yaml_node.IsDefined()) return;
    v_t<str_t> node_names;
    node_names.reserve(m_children.size());
    for (auto child: m_children) node_names.push_back(child->m_path.m_list.back());
    for (auto& content_name: child_names_in_file()) {
        auto found = std::find(node_names.cbegin(), node_names.cend(), content_name) != node_names.cend();
        REQUIRE_TRUE(found, log::format("file node \"{}\" is unrecognized", (m_path+content_name).m_string));
    }
}

v_t<std::pair<str_t, str_t>> conf_components::Node::help_pairs() const {
    return {};
}

bool conf_components::Node::make_in_file() const {
    if (!m_parent) return true;
    if (!m_parent->m_in_file) return false;
    try {
        m_parent->m_yaml_node[m_path.m_list.back()].Type();
        return true;
    }
    catch (const YAML::InvalidNode&){}
    return false;
}

conf_components::Node::Node(conf_components::Node* parent, const str_t& name, str_t desc) :
        m_path(parent->m_path+name), m_parent(parent), m_in_file(make_in_file()),
        m_yaml_node(m_in_file ? parent->m_yaml_node[name] : YAML::Node()),
        m_desc(std::move(desc)){
    if (parent) parent->m_children.push_back(this);
    if (m_in_file)
        REQUIRE_FALSE(m_yaml_node.IsNull(), "there should be no null values in YAML file");
}

conf_components::Node::Node(const YAML::Node& yaml_node, str_t desc) :
        m_path({}), m_parent(nullptr), m_in_file(true),
        m_yaml_node(yaml_node), m_desc(std::move(desc)){}

void conf_components::Node::validate() {
    for (auto child : m_children) child->validate();
    validate_file_contents();
    validate_node_contents();
}

void conf_components::Node::print_help(bool emph_first, size_t ilevel) const {
    auto pairs = help_pairs();
    if (!pairs.empty() && emph_first) pairs.front().second = log::bold_format(pairs.front().second);
    for (auto& pair: pairs) {
        std::cout << log::format(
                "{}{:<16}| {}\n",std::string(ilevel*2, ' '), pair.first, pair.second);
    }
    std::cout << '\n';
    for (auto child: m_children) child->print_help(emph_first, ilevel + 1);
}

void conf_components::Node::log() const {
    for (auto child: m_children) child->log();
}

str_t conf_components::enable_policy_string(conf_components::EnablePolicy enable_policy) {
    switch (enable_policy) {
        case Implicit: return "implicit";
        case Required: return "required";
        case Explicit: return "explicit";
    }
    return {};
}

bool conf_components::Group::make_enabled() const {
    /*
     * a Group without a parent is the root Node (Document) this is always enabled
     */
    if (!m_parent) return true;
    /*
     * if absent from file, only enable if Required or Implicit
     */
    if (!m_in_file) return m_enable_policy!=Explicit;
    /*
     * m_parent is a Node, convert it to Group
     */
    auto parent_group = dynamic_cast<const Group*>(m_parent);
    REQUIRE_TRUE(parent_group, "parent must be a group");
    /*
     * check that a non-Explicit does not have an Explicit parent
     */
    if (m_enable_policy!=Explicit)
        REQUIRE_FALSE(parent_group->m_enable_policy==Explicit,
                      "group cannot be implicitly enabled or required if its parent is explicitly enabled");
    /*
     * if the Node is in the file, and has a non-zero number of key-value pairs defined within it, it is enabled
     */
    if (m_yaml_node.IsMap()) return true;
    /*
     * otherwise, the group should be scalar, and boolean
     */
    REQUIRE_TRUE(m_yaml_node.IsScalar() && m_is_bool, "non-map type sections must be scalar boolean");
    /*
     * deal with the case that the group has been set to a scalar bool (true)
     */
    if (parse_as<bool>()) {
        if (m_enable_policy!=Explicit)
            log::warn("explicitly enabling {} group {} by boolean value is redundant",
                      m_enable_policy==Implicit ? "implicitly enabled" : "required", m_path.m_string);
        return true;
    }
    /*
     * we have group set to scalar bool (false), check that this is allowed by the policy before signaling disablement
     */
    REQUIRE_FALSE(m_enable_policy == Required, "required group cannot be disabled");
    return false;
}

conf_components::Group::Group(conf_components::Group* parent, const str_t& name, str_t desc,
                              conf_components::EnablePolicy enable_policy) :
        Node(parent, name, std::move(desc)), m_enable_policy(enable_policy),
        m_is_bool(m_in_file && m_yaml_node.IsScalar()), m_enabled(make_enabled()) {
    if (m_is_bool) REQUIRE_TRUE(parseable_as<bool>(),
                                "If group is scalar rather than map typed, it must be boolean");
}

conf_components::Group::Group(const YAML::Node& root, str_t desc) :
        Node(root, std::move(desc)), m_enable_policy(Required), m_is_bool(false), m_enabled(true){}

conf_components::Section::Section(conf_components::Group* parent, const str_t& name, str_t desc,
                                  conf_components::EnablePolicy enable_policy) :
        Group(parent, name, std::move(desc), enable_policy){}

v_t<std::pair<str_t, str_t>> conf_components::Section::help_pairs() const {
    return {
            {"Section",     m_path.m_string},
            {"Enabled",     enable_policy_string(m_enable_policy)},
            {"Description", m_desc}
    };
}

void conf_components::Section::log() const {
    log::info("{} ({})", m_path.m_string, m_enabled ? "enabled" : "disabled");
    if (m_enabled) Node::log();
}

bool conf_components::Document::contains_tabs(const str_t& contents) {
    std::regex r("(\\t)");
    auto begin = std::sregex_iterator(contents.cbegin(), contents.cend(), r);
    auto end = std::sregex_iterator();
    return std::distance(begin, end);
}

YAML::Node conf_components::Document::load(const str_t& fname) {
    /*
     * no file given: fills every Param with default values
     */
    if (fname.empty()) return {};
    str_t contents;
    /*
     * read in YAML file contents on the root node then then bcast to others
     */
    if (mpi::i_am_root()) contents = FileReader::to_string(fname);
    mpi::bcast(contents);
    REQUIRE_FALSE_ALL(contains_tabs(contents),
                      "tab characters are forbidden by the YAML standard: please use spaces");
    try {
        return YAML::Load(contents);
    }
    catch (const YAML::ParserException& ex) {
        ABORT(log::format("YAML syntax error in file {}, line {}, column {}",
                          fname, ex.mark.line, ex.mark.pos));
    }
    return {};
}

conf_components::Document::Document(const str_t& fname, str_t desc) : Group(load(fname), std::move(desc)) {}

conf_components::ParamBase::ParamBase(conf_components::Group* parent, const str_t& name, str_t desc, str_t dim_type,
                                      str_t default_value) :
        Node(parent, name, std::move(desc)),
        m_dim_type_str(std::move(dim_type)), m_default_value(std::move(default_value)){}

v_t<std::pair<str_t, str_t>> conf_components::ParamBase::help_pairs() const {
    return {
            {"Parameter",    m_path.m_string},
            {"Type",         m_dim_type_str},
            {"Default value",m_default_value},
            {"Description",  m_desc},
    };
}
