//
// Created by Robert J. Anderson on 25/06/2021.
//

#ifndef M7_YAMLWRAPPER_H
#define M7_YAMLWRAPPER_H

#include <yaml-cpp/yaml.h>
#include <yaml-cpp/node/node.h>
#include <yaml-cpp/parser.h>

#include <M7_lib/parallel/MPIAssert.h>

#include "Logging.h"
#include "FileReader.h"

namespace yaml {
    struct Path {
        std::list<std::string> m_name_list;

        Path(std::list<std::string> name_list);

        Path(std::vector<std::string> name_list);

        Path(std::string name);

        Path(const Path &other);

        std::string to_string() const;

        Path operator+(const std::string &name) const;

        Path up() const;

        size_t depth() const;
    };

    struct File {
        const std::string m_fname;
        YAML::Node m_root;

        File(const std::string &fname);

        YAML::Node get(Path path) const;

        YAML::Node get(std::list<std::string> path) const;

        YAML::Node get(std::string path) const;

        bool exists(Path path) const;

        bool exists(std::list<std::string> path) const;

        bool exists(std::string path) const;

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


#endif //M7_YAMLWRAPPER_H
