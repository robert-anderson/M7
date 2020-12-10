//
// Created by rja on 09/12/2020.
//

#ifndef M7_SERIALIZATION_H
#define M7_SERIALIZATION_H

#include "fstream"
#include <src/core/hash/Hashing.h>
#include "src/core/parallel/MPIWrapper.h"


struct Serializer {
    std::unique_ptr<std::ifstream> m_infile = nullptr;
    std::unique_ptr<std::ofstream> m_outfile = nullptr;
    std::vector<defs::data_t> m_buffer;
    defs::data_t* m_it;

    void reset_buffer_ptr() {
        m_it = m_buffer.data();
    }

    Serializer(size_t buffer_dsize): m_buffer(buffer_dsize){
        reset_buffer_ptr();
    }

    void set_input(std::string fname) {
        m_infile = std::unique_ptr<std::ifstream>(new std::ifstream(fname, std::ios::in | std::ios::binary));
    }

    void set_output(std::string fname) {
        m_outfile = std::unique_ptr<std::ofstream>(new std::ofstream(fname, std::ios::out | std::ios::binary));
    }

    size_t dwords_remaining() const {
        return (m_buffer.data()+m_buffer.size())-m_it;
    }

    void write() {
        if (!m_outfile) mpi::stop_all("Buffer to be flushed but no ofstream to write to");
        m_outfile->write((char*)m_buffer.data(), (m_it-m_buffer.data())*defs::nbyte_data);
        reset_buffer_ptr();
    }

    void read() {
        if (!m_infile) mpi::stop_all("Buffer to be filled but no ifstream to read from");
        m_infile->read((char*)m_buffer.data(), (m_it-m_buffer.data())*defs::nbyte_data);
        reset_buffer_ptr();
    }

    void save(char* src, size_t dsize){
        while (dsize){
            size_t dwords_to_copy = std::min(dsize, dwords_remaining());
            std::memcpy((char*)m_it, src, dwords_to_copy*defs::nbyte_data);
            m_it+=dwords_to_copy;
            dsize-=dwords_to_copy;
            if (!dwords_remaining()) write();
        }
    }

    template<typename T>
    void save(T v){
        save((char*)&v, mpi_dsizes[mpi_type_ind<T>()]);
    }

    void load(char* dst, size_t dsize){
        while (dsize){
            size_t dwords_to_copy = std::min(dsize, dwords_remaining());
            std::memcpy(dst, m_it, dwords_to_copy*defs::nbyte_data);
            m_it+=dwords_to_copy;
            dsize-=dwords_to_copy;
            if (!dwords_remaining()) read();
        }
    }

    template<typename T>
    T load(){
        T v;
        load((char*)&v, mpi_dsizes[mpi_type_ind<T>()]);
        return v;
    }

};


struct Serializable {

    struct Element {
        char* m_item_ptr;
        size_t m_type_ind;
        bool m_is_vector;

    private:
        template<typename T>
        const std::vector<T>* convert_vector() const{
            return reinterpret_cast<const std::vector<T>*>(m_item_ptr);
        }

        template<typename T>
        std::vector<T>* convert_vector() {
            return reinterpret_cast<std::vector<T>*>(m_item_ptr);
        }

        bool is_primitive() const {
            return m_type_ind!=~0ul;
        }

        char* vector_data() const {
            ASSERT(m_is_vector);
            switch (m_type_ind) {
                case 0: return (char*)convert_vector<char>()->data();
                case 1: return (char*)convert_vector<short int>()->data();
                case 2: return (char*)convert_vector<int>()->data();
                case 3: return (char*)convert_vector<long int>()->data();
                case 4: return (char*)convert_vector<long long int>()->data();
                case 5: return (char*)convert_vector<unsigned char>()->data();
                case 6: return (char*)convert_vector<unsigned short int>()->data();
                case 7: return (char*)convert_vector<unsigned int>()->data();
                case 8: return (char*)convert_vector<unsigned long int>()->data();
                case 9: return (char*)convert_vector<unsigned long long int>()->data();
                case 10: return (char*)convert_vector<float>()->data();
                case 11: return (char*)convert_vector<double>()->data();
                case 12: return (char*)convert_vector<long double>()->data();
                case 13: return (char*)convert_vector<std::complex<float>>()->data();
                case 14: return (char*)convert_vector<std::complex<double>>()->data();
                case 15: return (char*)convert_vector<std::complex<long double>>()->data();
            }
            return nullptr;
        }

        size_t vector_size() const {
            ASSERT(m_is_vector);
            switch (m_type_ind) {
                case 0: return convert_vector<char>()->size();
                case 1: return convert_vector<short int>()->size();
                case 2: return convert_vector<int>()->size();
                case 3: return convert_vector<long int>()->size();
                case 4: return convert_vector<long long int>()->size();
                case 5: return convert_vector<unsigned char>()->size();
                case 6: return convert_vector<unsigned short int>()->size();
                case 7: return convert_vector<unsigned int>()->size();
                case 8: return convert_vector<unsigned long int>()->size();
                case 9: return convert_vector<unsigned long long int>()->size();
                case 10: return convert_vector<float>()->size();
                case 11: return convert_vector<double>()->size();
                case 12: return convert_vector<long double>()->size();
                case 13: return convert_vector<std::complex<float>>()->size();
                case 14: return convert_vector<std::complex<double>>()->size();
                case 15: return convert_vector<std::complex<long double>>()->size();
            }
            return ~0ul;
        }

        size_t vector_dsize() const {
            ASSERT(is_primitive());
            return integer_utils::divceil(mpi_sizes[m_type_ind]*vector_size(), defs::nbyte_data);
        }

        size_t vector_capacity() const {
            ASSERT(m_is_vector);
            switch (m_type_ind) {
                case 0: return convert_vector<char>()->capacity();
                case 1: return convert_vector<short int>()->capacity();
                case 2: return convert_vector<int>()->capacity();
                case 3: return convert_vector<long int>()->capacity();
                case 4: return convert_vector<long long int>()->capacity();
                case 5: return convert_vector<unsigned char>()->capacity();
                case 6: return convert_vector<unsigned short int>()->capacity();
                case 7: return convert_vector<unsigned int>()->capacity();
                case 8: return convert_vector<unsigned long int>()->capacity();
                case 9: return convert_vector<unsigned long long int>()->capacity();
                case 10: return convert_vector<float>()->capacity();
                case 11: return convert_vector<double>()->capacity();
                case 12: return convert_vector<long double>()->capacity();
                case 13: return convert_vector<std::complex<float>>()->capacity();
                case 14: return convert_vector<std::complex<double>>()->capacity();
                case 15: return convert_vector<std::complex<long double>>()->capacity();
            }
            return ~0ul;
        }

        void vector_reserve(size_t n) {
            ASSERT(m_is_vector);
            switch (m_type_ind) {
                case 0: convert_vector<char>()->reserve(n); return;
                case 1: convert_vector<short int>()->reserve(n); return;
                case 2: convert_vector<int>()->reserve(n); return;
                case 3: convert_vector<long int>()->reserve(n); return;
                case 4: convert_vector<long long int>()->reserve(n); return;
                case 5: convert_vector<unsigned char>()->reserve(n); return;
                case 6: convert_vector<unsigned short int>()->reserve(n); return;
                case 7: convert_vector<unsigned int>()->reserve(n); return;
                case 8: convert_vector<unsigned long int>()->reserve(n); return;
                case 9: convert_vector<unsigned long long int>()->reserve(n); return;
                case 10: convert_vector<float>()->reserve(n); return;
                case 11: convert_vector<double>()->reserve(n); return;
                case 12: convert_vector<long double>()->reserve(n); return;
                case 13: convert_vector<std::complex<float>>()->reserve(n); return;
                case 14: convert_vector<std::complex<double>>()->reserve(n); return;
                case 15: convert_vector<std::complex<long double>>()->reserve(n); return;
            }
        }

        void vector_resize(size_t n) {
            ASSERT(m_is_vector);
            switch (m_type_ind) {
                case 0: convert_vector<char>()->resize(n); return;
                case 1: convert_vector<short int>()->resize(n); return;
                case 2: convert_vector<int>()->resize(n); return;
                case 3: convert_vector<long int>()->resize(n); return;
                case 4: convert_vector<long long int>()->resize(n); return;
                case 5: convert_vector<unsigned char>()->resize(n); return;
                case 6: convert_vector<unsigned short int>()->resize(n); return;
                case 7: convert_vector<unsigned int>()->resize(n); return;
                case 8: convert_vector<unsigned long int>()->resize(n); return;
                case 9: convert_vector<unsigned long long int>()->resize(n); return;
                case 10: convert_vector<float>()->resize(n); return;
                case 11: convert_vector<double>()->resize(n); return;
                case 12: convert_vector<long double>()->resize(n); return;
                case 13: convert_vector<std::complex<float>>()->resize(n); return;
                case 14: convert_vector<std::complex<double>>()->resize(n); return;
                case 15: convert_vector<std::complex<long double>>()->resize(n); return;
            }
        }

    public:
        size_t dsize() const {
            size_t res = 0ul;
            if(!m_is_vector) {
                if(!is_primitive()){
                    // scalar Serializable
                    res+=reinterpret_cast<Serializable*>(m_item_ptr)->dsize();
                } else {
                    // primative scalar
                    res+=mpi_dsizes[m_type_ind];
                }
            } else {
                // reserve 2 datawords for capacity and size
                res+=2;
                if(!is_primitive()){
                    // vector Serializable
                    auto ptr = reinterpret_cast<std::vector<Serializable>*>(m_item_ptr);
                    for (auto &it : *ptr) res+=it.dsize();
                } else {
                    // primative vector
                    res+=vector_dsize();
                }
            }
            return res;
        }

        size_t format_hash() const {
            defs::hash_t res = 0ul;
            if(!m_is_vector) {
                if(!is_primitive()){
                    // scalar Serializable
                    res = reinterpret_cast<Serializable*>(m_item_ptr)->format_hash();
                } else {
                    // primative scalar
                    res = mpi_sizes[m_type_ind];
                }
            } else {
                if(!is_primitive()){
                    // vector Serializable
                    auto ptr = reinterpret_cast<std::vector<Serializable>*>(m_item_ptr);
                    for (auto &it : *ptr) res^=it.format_hash();
                } else {
                    // primative vector
                    res = mpi_sizes[m_type_ind]+mpi_sizes.size();
                }
            }
            return res;
        }

        void save(Serializer& s) {
            if(!m_is_vector) {
                if(!is_primitive()){
                    // scalar Serializable
                    reinterpret_cast<Serializable*>(m_item_ptr)->save(s);
                } else {
                    // primative scalar
                    s.save(m_item_ptr, mpi_dsizes[m_type_ind]);
                }
            } else {
                if(!is_primitive()){
                    // vector Serializable
                    auto ptr = reinterpret_cast<std::vector<Serializable>*>(m_item_ptr);
                    s.save(ptr->capacity());
                    s.save(ptr->size());
                    for (auto &it : *ptr) {
                        it.save(s);
                    }
                } else {
                    // primative vector
                    s.save(vector_capacity());
                    s.save(vector_size());
                    s.save(vector_data(), vector_dsize());
                }
            }
        }


        void load(Serializer& s) {
            if(!m_is_vector) {
                if(!is_primitive()){
                    // scalar Serializable
                    reinterpret_cast<Serializable*>(m_item_ptr)->load(s);
                } else {
                    // primative scalar
                    s.load(m_item_ptr, mpi_dsizes[m_type_ind]);
                }
            } else {
                if(!is_primitive()){
                    // vector Serializable
                    auto ptr = reinterpret_cast<std::vector<Serializable>*>(m_item_ptr);
                    ptr->reserve(s.load<size_t>());
                    ptr->resize(s.load<size_t>());
                    for (auto &it : *ptr) {
                        it.load(s);
                    }
                } else {
                    vector_reserve(s.load<size_t>());
                    vector_resize(s.load<size_t>());
                    // primative vector
                    s.load(vector_data(), vector_dsize());
                }
            }
        }
    };

    std::vector<Element> m_elements;

    template<typename T>
    typename std::enable_if<mpi_supported<T>(), void>::type
    store(const T& data){
        auto type_ind = mpi_type_ind<T>();
        if (type_ind==~0ul) mpi::stop_all("Data is not serializable");
        m_elements.push_back({(char*)(&data), type_ind, false});
    }

    template<typename T>
    typename std::enable_if<mpi_supported<T>(), void>::type
    store(const std::vector<T>& data){
        auto type_ind = mpi_type_ind<T>();
        if (type_ind==~0ul) mpi::stop_all("Data is not serializable");
        m_elements.push_back({(char*)(&data), type_ind, true});
    }

    void store(const Serializable& data) {
        m_elements.push_back({(char*)(&data), ~0ul, false});
    }

    void store(const std::vector<Serializable>& data) {
        m_elements.push_back({(char*)(&data), ~0ul, true});
    }

    size_t dsize() const {
        size_t res = 0ul;
        for (auto element : m_elements) res +=element.dsize();
        return res;
    }

    /**
     * A checksum to verify that the data being read corresponds in format to the data being stored
     * @return
     */
    defs::hash_t format_hash() const {
        defs::hash_t res = 0;
        for (auto element : m_elements) {
            res ^= hashing::fnv_hash(element.format_hash());
            res <<=1;
        }
        return res;
    }

    void save(Serializer& s){
        for (auto element : m_elements) element.save(s);
    }

    void load(Serializer& s){
        for (auto element : m_elements) element.load(s);
    }

    void save_write(Serializer& s){
        s.reset_buffer_ptr();
        save(s);
        s.write();
    }

    void read_load(Serializer& s){
        s.read();
        load(s);
        s.reset_buffer_ptr();
    }


//    virtual char* serialize(char* ptr) {
//        for (auto pair : m_data){
//            if (pair.second==null_size)
//                ptr = reinterpret_cast<Serializable*>(pair.first)->serialize(ptr);
//            else {
//                std::memcpy(ptr, pair.first, pair.second);
//                ptr+=pair.second;
//            }
//        }
//        return ptr;
//    }
//
//    virtual const char* deserialize(const char* ptr){
//        for (auto pair : m_data){
//            if (pair.second==null_size)
//                ptr = reinterpret_cast<Serializable*>(pair.first)->deserialize(ptr);
//            else {
//                std::memcpy(pair.first, ptr, pair.second);
//                ptr+=pair.second;
//            }
//        }
//        return ptr;
//    }
};


#endif //M7_SERIALIZATION_H
