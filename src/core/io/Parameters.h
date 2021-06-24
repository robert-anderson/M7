//
// Created by rja on 23/06/2021.
//

#ifndef M7_ParamS_H
#define M7_ParamS_H

#include <src/core/parallel/MPIAssert.h>
#include "external/yaml-cpp/include/yaml-cpp/yaml.h"
#include "external/yaml-cpp/include/yaml-cpp/node/node.h"
#include "external/yaml-cpp/include/yaml-cpp/parser.h"
#include "Logging.h"
#include "FileReader.h"

namespace yaml {
    struct Path {
        std::list<std::string> m_path;

        Path(std::list<std::string> path) : m_path(std::move(path)) {}

        Path(std::vector<std::string> path) : m_path(path.cbegin(), path.cend()) {}

        Path(std::string path) : Path(string_utils::split(path, '.')) {}

        Path(const Path &other) : m_path(other.m_path) {}

        std::string to_string() const {
            auto it = m_path.cbegin();
            std::string out = *it;
            ++it;
            for (; it != m_path.cend(); ++it) out += "." + *it;
            return out;
        }

        Path operator+(const std::string &name) const {
            auto tmp = *this;
            tmp.m_path.push_back(name);
            return tmp;
        }

        Path up() const {
            auto tmp = *this;
            tmp.m_path.pop_back();
            return tmp;
        }

        size_t depth() const {
            return m_path.size();
        }
    };

    struct File {
        const std::string m_fname;
        YAML::Node m_root;

        File(const std::string &fname) : m_fname(fname) {
            std::string contents;
            /*
             * read in YAML file contents on the root node then then bcast to others
             */
            if (mpi::i_am_root()) contents = FileReader::to_string(fname);
            mpi::bcast(contents);
            try {
                m_root = YAML::Load(contents);
            }
            catch (YAML::ParserException ex) {
                ABORT(log::format("YAML syntax error in file {}", m_fname));
            }
        }

        YAML::Node get(Path path) const {
            auto node = m_root;
            for (const auto &it : path.m_path) node.reset(node[it]);
            return node;
        }

        YAML::Node get(std::list<std::string> path) const {
            return get(Path{path});
        }

        YAML::Node get(std::string path) const {
            return get(Path{path});
        }

        bool exists(Path path) const {
            return get(path).IsDefined();
        }

        bool exists(std::list<std::string> path) const {
            return exists(Path{path});
        }

        bool exists(std::string path) const {
            return exists(Path{path});
        }

        template<typename T>
        T get_as(const Path &path) const {
            return get(path).as<T>();
        }

        template<typename T>
        T get_as(std::list<std::string> path) const {
            return get(Path{path}).as<T>();
        }

        template<typename T>
        T get_as(std::string path) const {
            return get(Path{path}).as<T>();
        }
    };

}

namespace config {

    namespace opts {
        template<typename T>
        static std::string type_name() {
            if (std::is_same<T, size_t>::value) return "integer";
            if (std::is_same<T, double>::value) return "integer";
            if (std::is_same<T, std::string>::value) return "integer";
            if (std::is_same<T, std::string>::value) return "integer";
        }

        static std::string type_name(const double &) { return "float"; }

        static std::string type_name(const std::complex<double> &) { return "complex float"; }

        static std::string type_name(const std::string &) { return "string"; }

        template<typename T>
        static std::string type_name(const std::vector<T> &) { return "1D " + type_name<T>() + " string"; }
    }


    struct Node {
        const Node *m_parent;
        const yaml::Path m_path;
        const std::string m_description;
        std::list<const Node *> m_children;
        const std::string m_indent;

        Node(Node *parent, std::string name, std::string description) :
                m_parent(parent), m_path(parent->m_path + name),
                m_description(description), m_indent(2*(m_path.depth()-1), ' ') {
        }

        // only for ParamRoot
        Node(std::string description) :
                m_parent(nullptr), m_path(""), m_description(description), m_indent() {}


        virtual std::string help_string() const {
            return "";
        }

        virtual const yaml::File *get_file() const {
            REQUIRE_TRUE(m_parent, "Non-root Nodes must have a parent");
            return m_parent->get_file();
        }

        virtual std::string invalid_file_key() const {
            return "";
        }

        const std::string& name() const {
            return m_path.m_path.back();
        }
    };

    struct Group : Node {
    private:
        std::set<std::string> make_file_keys() const {
            std::set<std::string> file_keys;
            auto yf = get_file();
            if (!yf) return file_keys;
            auto yaml_node = yf->get(m_path);
            for (auto it = yaml_node.begin(); it != yaml_node.end(); ++it)
                file_keys.insert(YAML::Dump(it->first));
            return file_keys;
        }

        std::set<std::string> make_child_keys() const {
            std::set<std::string> child_keys;
            for (auto child: m_children) {
                auto child_key = child->name();
                REQUIRE_FALSE(child_keys.count(child_key), "Shouldn't have two child Nodes with the same name");
                child_keys.insert(child_key);
            }
            return child_keys;
        }

    public:
        Group(Group *parent, std::string name, std::string description) :
                Node(parent, name, description) {
            REQUIRE_TRUE(m_parent, "Non-root config::Nodes must have a parent");
            parent->add_child(this);
        }

        Group(std::string description) : Node(description) {}

        void add_child(const Node *child) {
            m_children.push_back(child);
        }

        std::string invalid_file_key() const override {
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
    };


    struct Section : Group {
        const bool m_is_specified;
    private:
        bool make_exists() const {
            auto file = get_file();
            if (file) return file->exists(m_path);
            return false;
        }

    public:
        Section(Group *parent, std::string name, std::string description) :
                Group(parent, name, description), m_is_specified(make_exists()) {}

        operator bool() const {
            return m_is_specified;
        }

        std::string help_string() const override {
            std::string str;
            str.append(log::format("{}Section:     {}\n", m_indent, m_path.to_string()));
            str.append(log::format("{}Description: {}\n\n", m_indent, m_description));
            for (auto child: m_children) str.append(child->help_string()+"\n");
            return str;
        }
    };

    struct Document : Group {
        const std::string m_name;
        const yaml::File *m_file;

        Document(const yaml::File *file, std::string name, std::string description) :
                Group(description), m_name(name), m_file(file) {}

        const yaml::File *get_file() const override {
            return m_file;
        }

        std::string help_string() const override {
            REQUIRE_FALSE(m_file, "Help string should only be generated when it has not been filled by a YAML file");
            std::string str = string_utils::boxed(log::format("Configuration document \"{}\"", m_name));
            str.append(log::format("Description: {}\n\n", m_description));
            for (auto child: m_children) str.append(child->help_string());
            return str;
        }
    };

    struct ParamBase : Node {
        const std::string m_v_default_str;
        const std::string m_dim_str;
        const std::string m_type_str;

        ParamBase(Group *parent, std::string name, std::string description, std::string v_default_str,
                  std::string dim_str, std::string type_str) :
                Node(parent, name, description), m_v_default_str(v_default_str),
                m_dim_str(dim_str), m_type_str(type_str) {
            REQUIRE_TRUE(m_parent, "Non-root config::Nodes must have a parent");
            parent->add_child(this);
        }

        std::string help_string() const override {
            std::string str;
            str.append(log::format("{}Parameter:     {}\n", m_indent, name()));
            str.append(log::format("{}Type:          {} {}\n", m_indent, m_dim_str, m_type_str));
            str.append(log::format("{}Default value: {}\n", m_indent, m_v_default_str));
            str.append(log::format("{}Description:   {}\n", m_indent, m_description));
            return str;
        }
    };

    template<typename T>
    static std::string type_str(const T&){
        ABORT(log::format("Unsupported type for a configuration parameter: {}",
                          log::get_demangled_symbol(typeid(T).name())));
        return "";
    }

    static std::string type_str(const size_t& v){
        return "integer";
    }

    static std::string type_str(const double& v){
        return "float";
    }

    static std::string type_str(const std::string& v){
        return "string";
    }

    template<typename T>
    static std::pair<std::string, std::string> dim_str(const T&){
        return "scalar";
    }

    template<typename T>
    static std::string dim_str(const std::vector<T>&){
        return "1D";
    }

    template<typename T>
    static std::string dim_str(const std::vector<std::vector<T>>&){
        return "2D";
    }

    template<typename T>
    static std::string dim_str(const std::vector<std::vector<std::vector<T>>>&){
        return "3D";
    }

    template<typename T>
    class Param : ParamBase {
        T m_v;

    public:
        Param(Group *parent, std::string name, std::string description, const T &v_default = {}):
                ParamBase(parent, name, description, utils::to_string(v_default),
                          dim_str(v_default), type_str(v_default)){
            auto file = parent->get_file();
            if (file) {
                try {
                    m_v = file->get_as<T>(m_path);
                }
                catch (YAML::BadConversion ex) {
                    ABORT(log::format("failed reading value {} from line {} of YAML config file",
                                      m_path.to_string(), ex.mark.line));
                }
            } else {
                m_v = v_default;
            }
        }

        const T &get() const {
            return m_v;
        }

        operator const T &() const {
            return m_v;
        }
    };
}


#endif //M7_ParamS_H
