//
// Created by rja on 21/10/2020.
//

#ifndef M7_FIELDSPECIFIER_H
#define M7_FIELDSPECIFIER_H

#include <cstddef>
#include <string>
#include <cstring>
#include <src/defs.h>
#include <map>
#include <src/core/hash/Hashing.h>
#include <src/core/parallel/MPIWrapper.h>

struct FieldData {
    const size_t m_element_size;
    const std::type_info &m_type_info;
    std::map<std::string, std::string> m_details;
};

struct FieldSpecifier {
    FieldData m_data;

    const size_t& element_size() const;

    const std::type_info& type_info() const;

    const std::vector<char> m_null_buffer;

    struct View {
        const FieldSpecifier &m_spec;
        char *m_ptr;

        View(const FieldSpecifier& field, char* ptr);

        View(const View &other);

        defs::data_t *dptr(const size_t &i) const;

        const size_t& element_size() const;

        void zero();

        int compare(const View& other) const;

        bool operator==(const View& other) const;

        bool operator!=(const View& other);

        bool is_zero() const;

        virtual std::string to_string() const = 0;

        void print() const;

        View &operator=(const View &other);

        void mpi_bcast(size_t iroot);

    };

    static defs::hash_t hash(const View& view);

    bool comparable_with(const FieldSpecifier& other) const;

    FieldSpecifier(size_t element_size, const std::type_info &type_info);

    virtual std::string element_string(char* ptr) const = 0;

};


#endif //M7_FIELDSPECIFIER_H
