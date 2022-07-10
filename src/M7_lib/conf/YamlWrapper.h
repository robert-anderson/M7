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

namespace conf_components {

    struct Path {
        const std::list<str_t> m_name_list;
        const str_t m_string;

        Path(std::list<str_t> name_list):
                m_name_list(name_list),
                m_string(string::join(v_t<str_t>(name_list.cbegin(), name_list.cend()), ".")){}

        Path operator+(const str_t &name) const {
            auto tmp = m_name_list;
            tmp.push_back(name);
            return {tmp};
        }

        Path up() const {
            auto tmp = m_name_list;
            tmp.pop_back();
            return tmp;
        }

        uint_t depth() const {
            return m_name_list.size();
        }
    };

    struct Node {
        const Path m_path;
        const bool m_impl_enabled;
    protected:
        const Node* m_parent;
        const YAML::Node m_yaml_node;
        const str_t m_desc;
        std::vector<Node*> m_children;

        v_t<str_t> child_names() const {
            if (!m_yaml_node.IsDefined()) return {};
            if (m_yaml_node.IsSequence()) return {};
            v_t<str_t> out;
            out.reserve(m_yaml_node.size());
            for (auto it=m_yaml_node.begin(); it!=m_yaml_node.end(); ++it) {
                out.push_back(it->first.Scalar());
            }
            return out;
        }

        template<typename T>
        bool parseable_as() const {
            try { m_yaml_node.as<T>(); }
            catch (const YAML::BadConversion &ex) { return false; }
            return true;
        }

        template<typename T>
        T parse_as() const {
            REQUIRE_TRUE(parseable_as<T>(), "value in file is not parseable as the required type");
            return m_yaml_node.as<T>();
        }

        bool null_in_file() const {
            if (!m_yaml_node.IsDefined() || m_yaml_node.size()) return false;
            return parse_as<str_t>() == "null";
        }

        Node(Path path, Node* parent, const YAML::Node& yaml_node, str_t desc, bool impl_enabled):
                m_path(std::move(path)), m_impl_enabled(impl_enabled),
                m_parent(parent), m_yaml_node(yaml_node), m_desc(std::move(desc)){
            if (parent) parent->m_children.push_back(this);
        }

        Node(const YAML::Node& yaml_node, str_t desc):
            Node({{}}, nullptr, yaml_node, std::move(desc), true){}

        virtual bool has_contents() const = 0;

        virtual str_t valid_logic() {return {};}

        /**
         * children may be defined in memory that are not present in the file contents, but file contents which do
         * not correspond to children in memory are invalid
         */
        str_t valid_contents() const {
            if (!m_yaml_node.IsDefined()) return {};
            v_t<str_t> memory_names;
            memory_names.reserve(m_children.size());
            for (auto child: m_children) memory_names.push_back(child->m_path.m_name_list.back());
            for (auto& content_name: child_names()) {
                if (std::find(memory_names.cbegin(), memory_names.cend(), content_name)==memory_names.cend())
                    return log::format("file node \"{}\" is unrecognized", (m_path+content_name).m_string);
            }
            return {};
        }

        virtual v_t<std::pair<str_t, str_t>> help_pairs() const {
            return {};
        }
    public:

        Node(Node* parent, const str_t& name, str_t desc, bool impl_enabled):
            Node(parent->m_path+name, parent, parent->m_yaml_node[name], std::move(desc), impl_enabled){
            REQUIRE_FALSE(null_in_file(), "there should be no null values in YAML file");
        }

        bool enabled() const {
            if (m_impl_enabled) return true;
            for (auto ptr = this; ptr; ptr=ptr->m_parent) {
                if (!ptr->m_impl_enabled && !has_contents()) return false;
            }
            return true;
        }

        str_t valid() {
            str_t str;
            str = valid_contents();
            if (!str.empty()) return str;
            str = valid_logic();
            if (!str.empty()) return str;
            for (auto child : m_children) {
                str = child->valid();
                if (!str.empty()) return str;
            }
            return str;
        }

        void require_valid() {
            auto str = valid();
            REQUIRE_TRUE(str.empty(), str);
        }

        void print_help(bool embolden_first = true, size_t ilevel=0) const {
            auto pairs = help_pairs();
            if (!pairs.empty() && embolden_first) pairs.front().second = log::bold_format(pairs.front().second);
            for (auto& pair: pairs) {
                std::cout << log::format(
                        "{}{:<16}| {}\n",std::string(ilevel*2, ' '), pair.first, pair.second);
            }
            std::cout << '\n';
            for (auto child: m_children) child->print_help(embolden_first, ilevel+1);
        }
    };

    struct Group : Node {

        Group(Node* parent, const str_t& name, str_t desc, bool impl_enabled=false):
                Node(parent, name, std::move(desc), impl_enabled){
            if (m_impl_enabled)
                REQUIRE_TRUE(m_parent->m_impl_enabled, "group cannot be implicitly enabled if its parent is not");
        }

    protected:

        Group(const YAML::Node& root, str_t desc): Node(root, std::move(desc)) {}

        /**
         * section must be defined in file and be either parseable as bool or contain children to be considered to have
         * a value in the file
         */
        bool has_contents() const override {
            if (!m_yaml_node.IsDefined()) return false;
            if (parseable_as<bool>()) return parse_as<bool>();
            return m_yaml_node.size();
        }

    };

    struct Section : Group {
        Section(Node* parent, const str_t& name, str_t desc, bool impl_enabled=false):
                Group(parent, name, std::move(desc), impl_enabled){}

        v_t<std::pair<str_t, str_t>> help_pairs() const override {
            return {
                    {"Section",     m_path.m_string},
                    {"Enable",      m_impl_enabled ? "implicit" : "explicit"},
                    {"Description", m_desc}
            };
        }
        // inherit implicit enabled from parent
        //Section(const Node* parent, const str_t& name): Section(parent, name, parent->m_impl_enabled){}
    };

    struct Document : Group {
    private:
        static bool contains_tabs(const str_t& contents){
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

    protected:
        bool has_contents() const override {
            return true;
        }
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
        ParamBase(Node* parent, const str_t& name, str_t desc, str_t dim_type, str_t default_value) :
            Node(parent, name, std::move(desc), parent->m_impl_enabled),
            m_dim_type_str(std::move(dim_type)), m_default_value(std::move(default_value)){}

        v_t<std::pair<str_t, str_t>> help_pairs() const override {
            return {
                    {"Parameter",    m_path.m_string},
                    {"Type",         m_dim_type_str},
                    {"Default value",m_dim_type_str},
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
                      convert::to_string(default_value)), m_default_value(default_value),
                      m_value(m_yaml_node.IsDefined() ? parse_as<T>() : m_default_value){}

    protected:
        bool has_contents() const override {
            return m_yaml_node.IsDefined() && parseable_as<T>();
        }

    public:
        operator const T&() const {
            return m_value;
        }

        Param& operator=(const T& v){
            m_value = v;
            return *this;
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
