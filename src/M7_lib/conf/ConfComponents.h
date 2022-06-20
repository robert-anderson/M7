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
        const std::string m_description;
        std::list<Node *> m_children;
        const std::string m_indent;

        Node(Node *parent, std::string name, std::string description);

        /**
         * only for ParamRoot
         */
        explicit Node(std::string description);

        virtual std::string help_string() const;

        virtual void log_value() const {}

        virtual void verify() {}

        virtual const yaml::File *get_file() const;

        virtual std::string invalid_file_key() const;

        const std::string &name() const;

        virtual bool enabled() const {
            return true;
        }

        bool parents_enabled() const;
    };

    struct Group : Node {
    private:
        std::set<std::string> make_file_keys() const;

        std::set<std::string> make_child_keys() const;

    public:
        Group(Group *parent, std::string name, std::string description);

        Group(std::string description);

        void add_child(Node *child);

        std::string invalid_file_key() const override;

        void verify() override {
            for (auto child: m_children) child->verify();
        }
    };


    struct Section : Group {
        Section(Group *parent, std::string name, std::string description);

        std::string help_string() const override;

        void log_value() const override {
            for (auto child: m_children) child->log_value();
        }
    };

    struct Document : Group {
        const std::string m_name;
        const yaml::File *m_file;

        Document(const yaml::File *file, std::string name, std::string description);

        const yaml::File *get_file() const override;

        std::string help_string() const override;

        void log_value() const override {
            log::info("Specified values for \"{}\"", m_name);
            for (auto child: m_children) child->log_value();
        }

        void verify() override;
    };

    struct ParamBase : Node {
        const std::string m_v_default_str;
        const std::string m_dim_type_str;

        ParamBase(Group *parent, std::string name, std::string description, std::string v_default_str,
                  std::string dim_type_str);

        std::string help_string() const override;
    };

    template<typename T=void>
    static std::string type_str() {
        ABORT(log::format("Unsupported type for a configuration parameter: {}",
                          log::get_demangled_symbol(typeid(T).name())));
        return "";
    }

    template<>
    std::string type_str<int>() { return "integer"; }

    template<>
    std::string type_str<long>() { return "long integer"; }

    template<>
    std::string type_str<size_t>() { return "unsigned integer"; }

    template<>
    std::string type_str<double>() { return "float"; }

    template<>
    std::string type_str<std::string>() { return "string"; }

    template<>
    std::string type_str<bool>() { return "boolean"; }

    template<typename T>
    static std::string dim_str(const T &) {
        return log::format("scalar {}", type_str<T>());
    }

    template<typename T>
    static std::string dim_str(const std::vector<T> &) {
        return log::format("1D {} array", type_str<T>());
    }

    template<typename T>
    static std::string dim_str(const std::vector<std::vector<T>> &) {
        return log::format("2D {} array", type_str<T>());
    }

    template<typename T>
    static std::string dim_str(const std::vector<std::vector<std::vector<T>>> &) {
        return log::format("3D {} array", type_str<T>());
    }

    template<typename T>
    class Param : public ParamBase {
        const T m_v_default;
        T m_v;

    public:
        Param(Group *parent, std::string name, const T &v_default, std::string description) :
                ParamBase(parent, name, description, utils::convert::to_string(v_default), dim_str(v_default)),
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
            log::info("{}: {}", m_yaml_path.to_string(), utils::convert::to_string(m_v));
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
