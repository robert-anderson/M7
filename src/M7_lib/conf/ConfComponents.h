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
        /**
         * pointer to parent Node (nullptr if document root)
         */
        const Node* m_parent;
        /**
         * true if the node exists in the file being parsed
         */
        const bool m_in_file;
    protected:
        /**
         * object of the yaml-cpp parser
         */
        const YAML::Node m_yaml_node;
        /**
         * text description displayed in help
         */
        const str_t m_desc;
        /**
         * all Nodes with this one as m_parent
         */
        std::vector<Node*> m_children;

    protected:
        /**
         * @return
         *  vector of all keys within this YAML Node
         */
        v_t<str_t> child_names_in_file() const;
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

        /**
         * check the loaded contents and perform any logical corrections to the tree (hence non-const)
         */
        virtual void validate_node_contents() {}

        /**
         * children may be defined in memory that are not present in the file contents, but file contents which do
         * not correspond to children in memory are invalid
         */
        void validate_file_contents() const;
        /**
         * @return
         *  help as a vector of key-value pairs
         */
        virtual v_t<std::pair<str_t, str_t>> help_pairs() const;

    private:
        /**
         * given that the m_parent and a m_path members are already set, determine whether this node is present in the
         * YAML input file
         * @return
         *  true if the node is present
         */
        bool make_in_file() const;

    protected:
        Node(Node* parent, const str_t& name, str_t desc);

        Node(const YAML::Node& yaml_node, str_t desc);

    public:
        /**
         * perform validation on all descendants, then on this Node. This means the modifications performed in
         * validate_node_contents occur in a bottom-up order
         */
        void validate();
        /**
         * recursively print the help to standard output
         * @param emph_first
         *  format the first value in bold
         * @param ilevel
         *  indentation level
         */
        void print_help(bool emph_first = true, size_t ilevel=0) const;
        /**
         * recursively log the state of the tree
         */
        virtual void log() const;
    };

    /**
     * enablement policy of a Group:
     *  - implicitly enabled: even if the Node does not appear in the input file, it is considered enabled
     *  - required: impliticly enabled, and cannot be disabled
     *  - explicitly enabled: considered disabled unless the Node appears in the input file
     */
    enum EnablePolicy {Implicit, Required, Explicit};

    str_t enable_policy_string(EnablePolicy enable_policy);

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
        /**
         * @return
         *  true if this Group is determined to be enabled
         */
        bool make_enabled() const;

    protected:
        Group(Group* parent, const str_t& name, str_t desc, EnablePolicy enable_policy=Required);

        Group(const YAML::Node& root, str_t desc);
    };

    struct Section : Group {
        Section(Group* parent, const str_t& name, str_t desc, EnablePolicy enable_policy=Required);

        v_t<std::pair<str_t, str_t>> help_pairs() const override;

        void log() const override;
    };


    /**
     * Group that only has Section children, and if enabled, only allows a certain number of those children to be
     * enabled
     */
    struct Selection : Group {
        enum Kind {AtLeast, NoMoreThan, Exactly};
        const Kind m_kind;
        const uint_t m_n;

        Selection(Group* parent, const str_t& name, str_t desc, Kind kind=NoMoreThan,
                  uint_t n=1ul, EnablePolicy enable_policy=Required);


    protected:
        void validate_node_contents() override;

        static str_t kind_string(Kind kind);

        v_t<std::pair<str_t, str_t>> help_pairs() const override;
    };

    struct Document : Group {
    private:
        static bool contains_tabs(const str_t& contents);

        static YAML::Node load(const str_t& fname);

    public:
        Document(const str_t& fname, str_t desc);
    };

    template<typename T=void>
    static str_t type_str() {
        ABORT(logging::format("Unsupported type for a configuration parameter: {}",
                          logging::get_demangled_symbol(typeid(T).name())));
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
        return logging::format("scalar {}", type_str<T>());
    }

    template<typename T>
    static str_t dim_str(const v_t<T> &) {
        return logging::format("1D {} array", type_str<T>());
    }

    template<typename T>
    static str_t dim_str(const v_t<v_t<T>> &) {
        return logging::format("2D {} array", type_str<T>());
    }

    template<typename T>
    static str_t dim_str(const v_t<v_t<v_t<T>>> &) {
        return logging::format("3D {} array", type_str<T>());
    }


    struct ParamBase : Node {
    private:
        const str_t m_dim_type_str, m_default_value;
    protected:
        ParamBase(Group* parent, const str_t& name, str_t desc, str_t dim_type, str_t default_value);

        v_t<std::pair<str_t, str_t>> help_pairs() const override;
    };

    template<typename T>
    struct Param : ParamBase {
        const T m_default_value;
        T m_value;

        Param(Group* parent, const str_t& name, T default_value, str_t desc):
            ParamBase(parent, name, std::move(desc), dim_str(default_value),
                      convert::to_string(default_value)), m_default_value(default_value){
            m_value = m_in_file ? parse_as<T>() : m_default_value;
        }

        operator const T&() const {
            return m_value;
        }

        Param& operator=(const T& v){
            m_value = v;
            return *this;
        }

        void log() const override {
            logging::info("{:<50}| {}", m_path.m_string, convert::to_string(m_value));
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
                out.push_back(logging::format("{} ({})", m_choices[i], m_choice_descriptions[i]));
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
            REQUIRE_TRUE(is_choices(v), logging::format("\"{}\" is not among the valid choices for param {}",
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
            tmp.emplace_back("Select one", ChoiceBase<T>::to_string());
            return tmp;
        }
    };

    template<typename T>
    class MultiChoice : public ChoiceBase<T>, Param<v_t<T>> {
        using ChoiceBase<T>::m_choices;
        using ChoiceBase<T>::m_choice_descriptions;
        using ChoiceBase<T>::are_choices;
        void require_are_choices(const v_t<T>& v) const {
            REQUIRE_TRUE(are_choices(v), logging::format("\"{}\" are not all among the valid choices for param {}",
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
            tmp.emplace_back("Select many", ChoiceBase<T>::to_string());
            return tmp;
        }
    };
}

#endif //M7_CONF_YAMLWRAPPER_H
