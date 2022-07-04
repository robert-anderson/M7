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
            return m_v;
        }

        void log_value() const override {
            log::info("{}: {}", m_yaml_path.to_string(), convert::to_string(m_v));
        }

        Param& operator=(const T& v){
            m_v = v;
            return *this;
        }
    };
}

template<typename T>
static std::ostream &operator<<(std::ostream &os, const conf_components::Param<T> &v) {
    os << v.get();
    return os;
}

#endif //M7_CONF_COMPONENTS_H
