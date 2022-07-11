//
// Created by rja on 10/07/22.
//

#ifndef M7_CONF_YAMLWRAPPER_H
#define M7_CONF_YAMLWRAPPER_H

#include <yaml-cpp/yaml.h>
#include <yaml-cpp/node/node.h>
#include <yaml-cpp/parser.h>

#include <M7_lib/parallel/MPIAssert.h>

#include <M7_lib/io/Logging.h>
#include <M7_lib/io/FileReader.h>
#include <regex>
#include <utility>

/**
 * class definitions related to input file parsing and configuring program behavior
 */
namespace conf_components {
    /**
     * sequence of string keys between document root and a specific Node in the YAML tree
     */
    struct Path {
        const std::list<str_t> m_list;
        const str_t m_string;

        Path(std::list<str_t> list);

        Path operator+(const str_t &name) const;

    };

    /**
     * basic object used to represent a Node in the tree structure of the YAML configuration document
     */
    struct Node {
        /**
         * absolute path from document root (YAML file Node) to this Node
         */
        const Path m_path;
    protected:
        /**
         * true if the node exists in the file being parsed
         */
        const bool m_in_file;
        /**
         * pointer to parent Node (nullptr if document root)
         */
        const Node* m_parent;
        /**
         * object of the yaml-cpp parser
         */
        const YAML::Node m_yaml_node;
        /**
         * text description displayed in help
         */
        const str_t m_desc;
        /**
         * all Nodes defined inside this one
         */
        std::vector<Node*> m_children;

    protected:
        /**
         * @return
         *  vector of all keys within this YAML Node
         */
        v_t<str_t> child_names_in_file() const {
            if (!m_yaml_node.IsDefined()) return {};
            if (m_yaml_node.IsSequence()) return {};
            v_t<str_t> out;
            out.reserve(m_yaml_node.size());
            for (auto it=m_yaml_node.begin(); it!=m_yaml_node.end(); ++it) {
                out.push_back(it->first.Scalar());
            }
            return out;
        }
        /**
         * @tparam T
         *  type to try and parse
         * @return
         *  true if yaml-cpp can parse the value as type T
         */
        template<typename T>
        bool parseable_as() const {
            try { m_yaml_node.as<T>(); }
            catch (const YAML::BadConversion &ex) { return false; }
            return true;
        }
        /**
         * @tparam T
         *  type to try and parse
         * @return
         *  parsed value
         */
        template<typename T>
        T parse_as() const {
            REQUIRE_TRUE(parseable_as<T>(), "value in file is not parseable as the required type");
            return m_yaml_node.as<T>();
        }

        bool in_file() const {
            try {
                m_yaml_node.Type();
                return true;
            }
            catch(const YAML::InvalidNode&){}
            return false;
        }
        /**
         * @return
         *  true if the file content associated with this Node has a null value e.g. "some_field: ~".
         *  such specifications are valid YAML but not accepted in this interface
         */
        bool null_in_file() const{return false;}

        /**
         * check the loaded contents and perform any logical corrections to the tree (hence non-const)
         */
        virtual void validate_node_contents() {}

        /**
         * children may be defined in memory that are not present in the file contents, but file contents which do
         * not correspond to children in memory are invalid
         */
        void validate_file_contents() const {
            if (!m_yaml_node.IsDefined()) return;
            v_t<str_t> node_names;
            node_names.reserve(m_children.size());
            for (auto child: m_children) node_names.push_back(child->m_path.m_list.back());
            for (auto& content_name: child_names_in_file()) {
                auto found = std::find(node_names.cbegin(), node_names.cend(), content_name) != node_names.cend();
                REQUIRE_TRUE(found, log::format("file node \"{}\" is unrecognized", (m_path+content_name).m_string));
            }
        }
        /**
         * @return
         *  help as a vector of key-value pairs
         */
        virtual v_t<std::pair<str_t, str_t>> help_pairs() const {
            return {};
        }

        Node(Node* parent, const str_t& name, str_t desc):
            m_path(parent->m_path+name), m_parent(parent),
            m_yaml_node(parent->m_yaml_node[name]), m_desc(std::move(desc)){
            if (parent) parent->m_children.push_back(this);
            REQUIRE_FALSE(null_in_file(), "there should be no null values in YAML file");
        }

        Node(YAML::Node yaml_node, str_t desc):
            m_path({}), m_parent(nullptr), m_yaml_node(std::move(yaml_node)), m_desc(std::move(desc)){}

    public:
        /**
         * perform validation on all descendants, then on this Node. This means the modifications performed in
         * validate_node_contents occur in a bottom-up order
         */
        void validate() {
            for (auto child : m_children) child->validate();
            validate_file_contents();
            validate_node_contents();
        }
        /**
         * recursively print the help to standard output
         * @param emph_first
         *  format the first value in bold
         * @param ilevel
         *  indentation level
         */
        void print_help(bool emph_first = true, size_t ilevel=0) const {
            auto pairs = help_pairs();
            if (!pairs.empty() && emph_first) pairs.front().second = log::bold_format(pairs.front().second);
            for (auto& pair: pairs) {
                std::cout << log::format(
                        "{}{:<16}| {}\n",std::string(ilevel*2, ' '), pair.first, pair.second);
            }
            std::cout << '\n';
            for (auto child: m_children) child->print_help(emph_first, ilevel + 1);
        }
        /**
         * recursively log the state of the tree
         */
        virtual void log() const {
            for (auto child: m_children) child->log();
        }
    };

    /**
     * enablement policy of a Group:
     *  - implicitly enabled: even if the Node does not appear in the input file, it is considered enabled
     *  - mandatory: impliticly enabled, and cannot be disabled
     *  - explicitly enabled: considered disabled unless the Node appears in the input file
     */
    enum EnablePolicy {Implicit, Mandatory, Explicit};

    str_t enable_policy_string(EnablePolicy enable_policy) {
        switch (enable_policy) {
            case Implicit: return "implicit";
            case Mandatory: return "mandatory";
            case Explicit: return "explicit";
        }
        return {};
    }

    struct Group : Node {
        const EnablePolicy m_enable_policy;
        /**
         * groups can be enabled and disabled by treating their key as that of a scalar boolean
         */
        const bool m_is_bool;
        /**
         * true if the state of the YAML file (given the enablement policy) implies that the functionality represented
         * by this Group
         */
        const bool m_enabled;

    private:
        bool make_enabled() const {

            if (!in_file()) return m_enable_policy!=Explicit;
            if (m_is_bool) REQUIRE_TRUE(parseable_as<bool>(),
                                        "If group is scalar rather than map typed, it must be boolean");

            if (!m_parent) return true;
            auto parent_group = dynamic_cast<const Group*>(m_parent);
            REQUIRE_TRUE(parent_group, "parent must be a group");
            if (m_enable_policy!=Explicit)
                REQUIRE_FALSE(parent_group->m_enable_policy==Explicit,
                              "group cannot be implicitly enabled or mandatory if its parent is explicitly enabled");
            if (m_yaml_node.size()) return true;
            REQUIRE_TRUE(m_is_bool, "null nodes are forbidden");
            if (parse_as<bool>()) {
                if (m_enable_policy!=Explicit)
                    log::warn("explicitly enabling {} group {} by boolean value is redundant",
                              m_enable_policy==Implicit ? "implicitly enabled" : "mandatory", m_path.m_string);
                return true;
            }
            REQUIRE_FALSE(m_enable_policy==Mandatory, "mandatory group cannot be disabled");
            return false;
        }

    protected:
        Group(Group* parent, const str_t& name, str_t desc, EnablePolicy enable_policy=Mandatory):
                Node(parent, name, std::move(desc)), m_enable_policy(enable_policy),
                m_is_bool(in_file() && m_yaml_node.IsScalar()), m_enabled(make_enabled()) {}

        Group(const YAML::Node& root, str_t desc) :
                Node(root, std::move(desc)), m_enable_policy(Mandatory), m_is_bool(false), m_enabled(true){}
    };

    struct Section : Group {
        Section(Group* parent, const str_t& name, str_t desc, EnablePolicy enable_policy=Mandatory) :
                Group(parent, name, std::move(desc), enable_policy){}

        v_t<std::pair<str_t, str_t>> help_pairs() const override {
            return {
                    {"Section",     m_path.m_string},
                    {"Enabled",     enable_policy_string(m_enable_policy)},
                    {"Description", m_desc}
            };
        }

        void log() const override {
            log::info("{} ({})", m_path.m_string, m_enabled ? "enabled" : "disabled");
            if (m_enabled) Node::log();
        }
    };

    struct Document : Group {
    private:
        static bool contains_tabs(const str_t& contents) {
            std::regex r("(\\t)");
            auto begin = std::sregex_iterator(contents.cbegin(), contents.cend(), r);
            auto end = std::sregex_iterator();
            return std::distance(begin, end);
        }

        static YAML::Node load(const str_t& fname) {
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

    public:
        Document(const str_t& fname, str_t desc): Group(load(fname), std::move(desc)) {}
    };

    template<typename T=void>
    static str_t type_str() {
        ABORT(log::format("Unsupported type for a configuration parameter: {}",
                          log::get_demangled_symbol(typeid(T).name())));
        return "";
    }

    template<>
    str_t type_str<int>() { return "integer"; }

    template<>
    str_t type_str<long>() { return "long integer"; }

    template<>
    str_t type_str<uint_t>() { return "unsigned integer"; }

    template<>
    str_t type_str<double>() { return "float"; }

    template<>
    str_t type_str<str_t>() { return "string"; }

    template<>
    str_t type_str<bool>() { return "boolean"; }

    template<typename T>
    static str_t dim_str(const T &) {
        return log::format("scalar {}", type_str<T>());
    }

    template<typename T>
    static str_t dim_str(const v_t<T> &) {
        return log::format("1D {} array", type_str<T>());
    }

    template<typename T>
    static str_t dim_str(const v_t<v_t<T>> &) {
        return log::format("2D {} array", type_str<T>());
    }

    template<typename T>
    static str_t dim_str(const v_t<v_t<v_t<T>>> &) {
        return log::format("3D {} array", type_str<T>());
    }


    struct ParamBase : Node {
    private:
        const str_t m_dim_type_str, m_default_value;
    protected:
        ParamBase(Group* parent, const str_t& name, str_t desc, str_t dim_type, str_t default_value) :
                Node(parent, name, std::move(desc)),
                m_dim_type_str(std::move(dim_type)), m_default_value(std::move(default_value)){}

        v_t<std::pair<str_t, str_t>> help_pairs() const override {
            return {
                    {"Parameter",    m_path.m_string},
                    {"Type",         m_dim_type_str},
                    {"Default value",m_default_value},
                    {"Description",  m_desc},
            };
        }
    };

    template<typename T>
    struct Param : ParamBase {
        const T m_default_value;
        T m_value;

        Param(Group* parent, const str_t& name, T default_value, str_t desc):
            ParamBase(parent, name, std::move(desc), dim_str(default_value),
                      convert::to_string(default_value)), m_default_value(default_value){
            m_value = in_file() ? parse_as<T>() : m_default_value;
        }

        operator const T&() const {
            return m_value;
        }

        Param& operator=(const T& v){
            m_value = v;
            return *this;
        }

        void log() const override {
            log::info("{:<40}| {}", m_path.m_string, convert::to_string(m_value));
        }
    };

    template<typename T>
    struct ChoiceBase {
        const v_t<T> m_choices;
        const v_t<str_t> m_choice_descriptions;

        ChoiceBase(v_t<T> choices, v_t<str_t> choice_descriptions = {}):
                m_choices(std::move(choices)), m_choice_descriptions(std::move(choice_descriptions)){
            if (!choice_descriptions.empty()) REQUIRE_EQ(m_choices.size(), m_choice_descriptions.size(),
                                                         "if descriptions are defined, each choice must have one");
        }
        ChoiceBase(const v_t<std::pair<T, str_t>>& choice_and_descs):
                ChoiceBase<T>(convert::split(choice_and_descs).first, convert::split(choice_and_descs).second){}

        str_t to_string() const {
            if (m_choice_descriptions.empty()) return convert::to_string(m_choices);
            v_t<str_t> out;
            out.reserve(m_choices.size());
            for (uint_t i=0ul; i<m_choices.size(); ++i)
                out.push_back(log::format("{} ({})", m_choices[i], m_choice_descriptions[i]));
            return "["+string::join(out, ", ")+"]";
        }
    protected:
        bool is_a_choice(const T& v) const {
            return std::find(m_choices.cbegin(), m_choices.cend(), v)!=m_choices.cend();
        }

        bool are_choices(const v_t<T>& v) const {
            return std::all_of(v.cbegin(), v.cend(), [this](const T& it){return is_a_choice(it);});
        }
    };

    template<typename T>
    class SingleChoice : public ChoiceBase<T>, Param<T> {
        using ChoiceBase<T>::m_choices;
        using ChoiceBase<T>::m_choice_descriptions;
        using ChoiceBase<T>::is_choices;
        void require_is_choice(const T& v) const {
            REQUIRE_TRUE(is_choices(v), log::format("\"{}\" is not among the valid choices for param {}",
                                                    convert::to_string(v), Node::m_path.m_string));
        }
        static T get_first(const v_t<T>& choices) {
            REQUIRE_FALSE(choices.empty(), "choices vector must be non-empty");
            return choices.front();
        }
        static T get_first(const v_t<std::pair<T, str_t>>& choice_and_descs) {
            REQUIRE_FALSE(choice_and_descs.empty(), "choices vector must be non-empty");
            return choice_and_descs.front().first;
        }
    public:
        using Param<T>::m_value;
        using Param<T>::operator const T&;
    private:
        void validate(const T &v_default) const {
            require_is_choice(v_default);
            require_is_choice(m_value);
        }

    private:
        SingleChoice(Group *parent, str_t name, v_t<T> choices, v_t<str_t> choice_descriptions,
                     const T &v_default, str_t description):
                ChoiceBase<T>(choices, choice_descriptions),
                Param<T>(parent, name, v_default, description){
            validate(v_default);
        }

    public:
        SingleChoice(Group *parent, str_t name, const v_t<std::pair<T, str_t>>& choice_and_descs,
                     const v_t<T> &v_default, str_t description):
                ChoiceBase<T>(choice_and_descs), Param<v_t<T>>(parent, name, v_default, description){
            validate(v_default);
        }

        SingleChoice(Group *parent, str_t name, v_t<T> choices, const v_t<T> &v_default, str_t description):
                ChoiceBase<T>(choices), Param<v_t<T>>(parent, name, v_default, description){
            validate(v_default);
        }
        /*
         * default value is taken to be the first choice
         */
        SingleChoice(Group *parent, str_t name, const v_t<std::pair<T, str_t>>& choice_and_descs, str_t description):
                SingleChoice(parent, name, choice_and_descs, get_first(choice_and_descs), description){}

        SingleChoice(Group *parent, str_t name, const v_t<T>& choices, str_t description):
                SingleChoice(parent, name, choices, get_first(choices), description){}

    private:

        v_t<std::pair<str_t, str_t>> help_pairs() const override {
            auto tmp = ParamBase::help_pairs();
            tmp.emplace_back("Select one from", ChoiceBase<T>::to_string());
            return tmp;
        }
    };

    template<typename T>
    class MultiChoice : public ChoiceBase<T>, Param<v_t<T>> {
        using ChoiceBase<T>::m_choices;
        using ChoiceBase<T>::m_choice_descriptions;
        using ChoiceBase<T>::are_choices;
        void require_are_choices(const v_t<T>& v) const {
            REQUIRE_TRUE(are_choices(v), log::format("\"{}\" are not all among the valid choices for param {}",
                                                     convert::to_string(v), Node::m_path.m_string));
        }
    public:
        using Param<v_t<T>>::m_value;
        using Param<v_t<T>>::operator const v_t<T>&;
    private:
        void validate(const v_t<T> &v_default, bool repetition) const {
            require_are_choices(v_default);
            require_are_choices(m_value);
            if (!repetition) {
                REQUIRE_EQ(std::set<T>(m_value.begin(), m_value.end()).size(),
                           m_value.size(), "repeated choices are not allowed");
            }
        }

    private:
        MultiChoice(Group *parent, str_t name, v_t<T> choices, v_t<str_t> choice_descriptions,
                    const T &v_default, str_t description, bool repetition):
                ChoiceBase<T>(choices, choice_descriptions),
                Param<v_t<T>>(parent, name, v_default, description){
            validate(v_default, repetition);
        }

    public:
        MultiChoice(Group *parent, str_t name, const v_t<std::pair<T, str_t>>& choice_and_descs,
                    const v_t<T> &v_default, str_t description, bool repetition):
                ChoiceBase<T>(choice_and_descs), Param<v_t<T>>(parent, name, v_default, description){
            validate(v_default, repetition);
        }

        MultiChoice(Group *parent, str_t name, v_t<T> choices, const v_t<T> &v_default,
                    str_t description, bool repetition):
                ChoiceBase<T>(choices), Param<v_t<T>>(parent, name, v_default, description){
            validate(v_default, repetition);
        }

    private:

        v_t<std::pair<str_t, str_t>> help_pairs() const override {
            auto tmp = ParamBase::help_pairs();
            tmp.emplace_back("Select many from", ChoiceBase<T>::to_string());
            return tmp;
        }
    };
}

#endif //M7_CONF_YAMLWRAPPER_H
