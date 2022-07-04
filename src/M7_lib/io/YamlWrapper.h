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
        std::list<str_t> m_name_list;

        Path(std::list<str_t> name_list);

        Path(strv_t name_list);

        Path(str_t name);

        Path(const Path &other);

        str_t to_string() const;

        Path operator+(const str_t &name) const;

        Path up() const;

        uint_t depth() const;
    };

    struct File {
        const str_t m_fname;
        YAML::Node m_root;

        File(const str_t &fname);

        YAML::Node get(Path path) const;

        YAML::Node get(std::list<str_t> path) const;

        YAML::Node get(str_t path) const;

        bool exists(Path path) const;

        bool exists(std::list<str_t> path) const;

        bool exists(str_t path) const;

        template<typename T>
        T get_as(const Path &path) const {
            return get(path).as<T>();
        }

        template<typename T>
        T get_as(std::list<str_t> path) const {
            return get(Path{path}).as<T>();
        }

        template<typename T>
        T get_as(str_t path) const {
            return get(Path{path}).as<T>();
        }
    };

}


#endif //M7_YAMLWRAPPER_H
