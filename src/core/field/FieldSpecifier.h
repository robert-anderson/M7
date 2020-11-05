//
// Created by rja on 21/10/2020.
//

#ifndef M7_FIELDSPECIFIER_H
#define M7_FIELDSPECIFIER_H

#include <cstddef>
#include <string>
#include <cstring>
#include <src/core/util/defs.h>
#include <map>
#include <src/core/hash/Hashing.h>

struct FieldData {
    const size_t m_element_size;
    const std::type_info &m_type_info;
    std::map<std::string, std::string> m_details;
};

struct FieldSpecifier {
    FieldData m_data;

    const size_t& element_size() const {
        return m_data.m_element_size;
    };
    const std::type_info& type_info() const {
        return m_data.m_type_info;
    };

    const std::vector<char> m_null_buffer;
    struct View {
        const FieldSpecifier &m_spec;
        char *m_ptr;

        View(const FieldSpecifier& field, char* ptr): m_spec(field), m_ptr(ptr){}

        defs::data_t *dptr(const size_t &i) const {
            ASSERT(i * defs::nbyte_data < m_spec.element_size());
            return ((defs::data_t *) m_ptr) + i;
        }

        const size_t& element_size() const {
            return m_spec.element_size();
        }

        void zero() {
            std::memset(m_ptr, 0, element_size());
        }

        int compare(const View& other) const {
            ASSERT(element_size()==other.element_size());
            return std::memcmp(m_ptr, other.m_ptr, element_size());
        }

        bool operator==(const View& other) const {
            return compare(other)==0;
        }

        bool operator!=(const View& other){
            return !(*this==other);
        }

        bool is_zero() const {
            return !std::memcmp(m_ptr, m_spec.m_null_buffer.data(), element_size());
        }

        virtual std::string to_string() const = 0;

        void print() const {
            std::cout << to_string() << std::endl;
        }

        View &operator=(const View &other) {
            ASSERT(m_spec.element_size() == other.m_spec.element_size());
            if (&other != this)
                std::memcpy(m_ptr, other.m_ptr, m_spec.element_size());
            return *this;
        }

    protected:

        View(const View &other) : m_spec(other.m_spec), m_ptr(other.m_ptr) {}

    };

    static defs::hash_t hash(const View& view) {
        return hashing::fnv_hash(view.m_ptr, view.element_size());
    }

    bool comparable_with(const FieldSpecifier& other) const {
        return (type_info()==other.type_info()) && (element_size()==other.element_size());
    }


//    template<typename ...Args>
//    defs::hash_t hash(const_raw_view_t first, Args... rest) const {
//        return hashing::fnv_hash(first.first, first.second)^hash(rest...);
//    }
//
//    template<typename ...Args>
//    defs::hash_t hash(const View& first, Args... rest) const {
//        ASSERT(comparable_with(first.m_spec));
//        return hash(convert_to_raw(first), rest...);
//    }

    FieldSpecifier(size_t element_size, const std::type_info &type_info);

    virtual std::string element_string(char* ptr) const = 0;

};


#endif //M7_FIELDSPECIFIER_H
