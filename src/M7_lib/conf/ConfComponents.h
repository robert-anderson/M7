//
// Created by Robert J. Anderson on 25/06/2021.
//

#ifndef M7_CONF_COMPONENTS_H
#define M7_CONF_COMPONENTS_H

#include <M7_lib/io/YamlWrapper.h>

namespace conf_components {

    struct Node {
        const Node *m_parent;
        const yaml::Path m_yaml_path;
        const str_t m_description;
        std::list<Node *> m_children;
        const str_t m_indent;

        Node(Node *parent, str_t name, str_t description);

        /**
         * only for ParamRoot
         */
        explicit Node(str_t description);

        virtual str_t help_string() const;

        virtual void log_value() const {}

        virtual void verify() {}

        virtual const yaml::File *get_file() const;

        virtual str_t invalid_file_key() const;

        const str_t &name() const;

        virtual bool enabled() const {
            return true;
        }

        bool parents_enabled() const;
    };

    struct Group : Node {
    private:
        std::set<str_t> make_file_keys() const;

        std::set<str_t> make_child_keys() const;

    public:
        Group(Group *parent, str_t name, str_t description);

        Group(str_t description);

        void add_child(Node *child);

        str_t invalid_file_key() const override;

        void verify() override {
            for (auto child: m_children) child->verify();
        }
    };


    struct Section : Group {
        Section(Group *parent, str_t name, str_t description);

        str_t help_string() const override;

        void log_value() const override {
            for (auto child: m_children) child->log_value();
        }
    };

    struct Document : Group {
        const str_t m_name;
        const yaml::File *m_file;

        Document(const yaml::File *file, str_t name, str_t description);

        const yaml::File *get_file() const override;

        str_t help_string() const override;

        void log_value() const override {
            log::info("Specified values for \"{}\"", m_name);
            for (auto child: m_children) child->log_value();
        }

        void verify() override;
    };

    struct ParamBase : Node {
        const str_t m_v_default_str;
        const str_t m_dim_type_str;

        ParamBase(Group *parent, str_t name, str_t description, str_t v_default_str,
                  str_t dim_type_str);

        str_t help_string() const override;
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

    template<typename T>
    class Param : public ParamBase {
        const T m_v_default;
    protected:
        T m_v;
    public:
        Param(Group *parent, str_t name, const T &v_default, str_t description) :
                ParamBase(parent, name, description, convert::to_string(v_default), dim_str(v_default)),
                m_v_default(v_default) {
            auto file = parent->get_file();
            if (file) {
                try {
                    if (file->exists(m_yaml_path)) m_v = file->get_as<T>(m_yaml_path);
                    else m_v = m_v_default;
                }
                catch (const YAML::BadConversion &ex) {
                    ABORT(log::format("failed reading value {} from line {} of YAML config file",
                                      m_yaml_path.to_string(), ex.mark.line));
                }
            } else {
                m_v = m_v_default;
            }
        }

        const T &get() const {
            return m_v;
        }

        operator const T &() const {
            return get();
        }

        void log_value() const override {
            log::info("{}: {}", m_yaml_path.to_string(), convert::to_string(m_v));
        }

        Param& operator=(const T& v){
            m_v = v;
            return *this;
        }
    };

#if 0
    template<typename T>
    class Choice : public Param<> {
    protected:
        const v_t<T> m_choices;
        const v_t<str_t> m_choice_descriptions;
        bool is_a_choice(const T& v) const {
            return std::find(m_choices.cbegin(), m_choices.cend(), v)!=m_choices.cend();
        }
        void require_is_a_choice(const T& v) const {
            REQUIRE_TRUE(is_a_choice(v), log::format("\"{}\" is not among the valid choices for param {}",
                                                     convert::to_string(v), Node::m_yaml_path.to_string()));
        }
        bool are_choices(const v_t<T>& v) {
            return std::all_of(v.cbegin(), v.cend(), [this](const T& it){return is_a_choice(it);});
        }
        void require_are_choices(const v_t<T>& v) {
            REQUIRE_TRUE(are_choices(v), log::format("\"{}\" are not all among the valid choices for param {}",
                                                     convert::to_string(v), Node::m_yaml_path.to_string()));
        }

    private:
        Choice(Group *parent, str_t name, v_t<T> choices, v_t<str_t> choice_descriptions,
               const T &v_default, str_t description): Param<T>(parent, name, v_default, description),
                    m_choices(std::move(choices)),
                    m_choice_descriptions(std::move(choice_descriptions)){}

    protected:
        Choice(Group *parent, str_t name, const v_t<std::pair<T, str_t>>& choice_and_descs,
               const T &v_default, str_t description):
               Choice(parent, name, convert::split(choice_and_descs).first,
                      convert::split(choice_and_descs).second, v_default, description){}

        Choice(Group *parent, str_t name, v_t<T> choices, const T &v_default, str_t description):
                Choice(parent, name, choices, {}, v_default, description){}

        str_t to_string() const {
            if (m_choice_descriptions.empty()) return convert::to_string(m_choices);
            v_t<str_t> out;
            out.reserve(m_choices.size());
            for (uint_t i=0ul; i<m_choices.size(); ++i)
                out.push_back(log::format("{} ({})", m_choices[i], m_choice_descriptions[i]));
            return convert::to_string(out);
        }
    };

    template<typename T>
    class SingleChoice : public Choice<T> {
        static T get_first(const v_t<T>& choices) {
            REQUIRE_FALSE(choices.empty(), "choices vector must be non-empty");
            return choices.front();
        }
        static T get_first(const v_t<std::pair<T, str_t>>& choice_and_descs) {
            REQUIRE_FALSE(choice_and_descs.empty(), "choices vector must be non-empty");
            return choice_and_descs.front().first;
        }
    public:
        SingleChoice(Group *parent, str_t name, const v_t<std::pair<T, str_t>>& choice_and_descs,
                     const T &v_default, str_t description): Choice<T>(parent, name, choice_and_descs, v_default, description){
            Choice<T>::require_is_a_choice(v_default);
            Choice<T>::require_is_a_choice(Param<T>::get());
        }

        SingleChoice(Group *parent, str_t name, const v_t<T>& choices, const T &v_default, str_t description):
                Choice<T>(parent, name, choices, v_default, description){
            Choice<T>::require_is_a_choice(v_default);
            Choice<T>::require_is_a_choice(Param<T>::get());
        }
        /*
         * default value is taken to be the first choice
         */
        SingleChoice(Group *parent, str_t name, const v_t<std::pair<T, str_t>>& choice_and_descs, str_t description):
                Choice<T>(parent, name, choice_and_descs, get_first(choice_and_descs), description){}

        SingleChoice(Group *parent, str_t name, const v_t<T>& choices, str_t description):
                Choice<T>(parent, name, choices, get_first(choices), description){}

        str_t help_string() const override {
            auto str = ParamBase::help_string();
            str.append(log::format("{}Select one from:  {}\n", ParamBase::m_indent, Choice<T>::to_string()));
            return str;
        }

        using Param<T>::get;
        using Param<T>::operator const T &;
    };

    template<typename T>
    class MultiChoice : public Choice<v_t<T>> {
    public:
        MultiChoice(Group *parent, str_t name, const v_t<std::pair<T, str_t>>& choice_and_descs,
                     const v_t<T>& v_default, str_t description):
            Choice<v_t<T>>(parent, name, choice_and_descs, v_default, description){
            Choice<T>::require_are_choices(v_default);
            Choice<T>::require_are_choices(Param<T>::get());
        }

        MultiChoice(Group *parent, str_t name, const v_t<T>& choices, const v_t<T>& v_default, str_t description):
            Choice<v_t<T>>(parent, name, choices, v_default, description){
            Choice<T>::require_are_choices(v_default);
            Choice<T>::require_are_choices(Param<T>::get());
        }

        using Param<T>::m_v;
        bool contains(const T& v) const {
            return std::find(m_v.cbegin(), m_v.cend(), v)!=m_v.cend();
        }

        str_t help_string() const override {
            auto str = ParamBase::help_string();
            str.append(log::format("{}Select many from: {}\n", ParamBase::m_indent, Choice<v_t<T>>::to_string()));
            return str;
        }
        using Param<T>::get;
        using Param<T>::operator const T &;
    };
#endif

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
    class MultiChoice : public ChoiceBase<T>, Param<v_t<T>> {
        using ChoiceBase<T>::m_choices;
        using ChoiceBase<T>::m_choice_descriptions;
        using ChoiceBase<T>::are_choices;
        void require_are_choices(const v_t<T>& v) const {
            REQUIRE_TRUE(are_choices(v), log::format("\"{}\" are not all among the valid choices for param {}",
                                                     convert::to_string(v), Node::m_yaml_path.to_string()));
        }
    public:
        using Param<v_t<T>>::get;
        using Param<v_t<T>>::operator const v_t<T>&;
    private:
        void validate(const v_t<T> &v_default, bool repetition) const {
            require_are_choices(v_default);
            require_are_choices(get());
            if (!repetition) {
                REQUIRE_EQ(std::set<T>(get().begin(), get().end()).size(), get().size(), "repeated choices are not allowed");
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
        str_t help_string() const override {
            auto str = ParamBase::help_string();
            str.append(log::format("{}Select many from: {}\n", ParamBase::m_indent, ChoiceBase<T>::to_string()));
            return str;
        }
    };
}

template<typename T>
static std::ostream &operator<<(std::ostream &os, const conf_components::Param<T> &v) {
    os << v.get();
    return os;
}

#endif //M7_CONF_COMPONENTS_H
