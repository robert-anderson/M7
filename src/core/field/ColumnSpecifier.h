//
// Created by rja on 21/10/2020.
//

#ifndef M7_COLUMNSPECIFIER_H
#define M7_COLUMNSPECIFIER_H

#include <cstddef>
#include <string>
#include <cstring>
#include <src/defs.h>
#include <map>
#include <src/core/hash/Hashing.h>
#include <src/core/parallel/MPIWrapper.h>

struct ColumnData {
    const size_t m_element_size;
    const std::type_info &m_type_info;
    std::map<std::string, std::string> m_details;
};

struct ColumnSpecifier {
    ColumnData m_data;

    struct View {
        const ColumnSpecifier &m_spec;
        char *m_ptr;
        View(const ColumnSpecifier& spec, char* ptr);
        View(const View &other);
        const size_t& element_size() const;
        int compare(const View& other) const;
        bool operator==(const View& other) const;
        bool operator!=(const View& other);
        bool is_zero() const;
        virtual std::string to_string() const = 0;
        void print() const;
        const defs::data_t *cdptr(const size_t &i) const;
        defs::data_t *dptr(const size_t &i) const;
        void zero();
        void mpi_bcast(size_t iroot);
        View &operator=(const View &other);
    };

    ColumnSpecifier(size_t element_size, const std::type_info &type_info);

    static defs::hash_t hash(const View& view);

    bool comparable_with(const ColumnSpecifier& other) const;

    virtual std::string element_string(char* ptr) const = 0;

    const size_t& element_size() const;

    const std::type_info& type_info() const;

    const std::vector<char> m_null_buffer;
};


#endif //M7_COLUMNSPECIFIER_H
