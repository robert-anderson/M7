//
// Created by rja on 23/06/22.
//

#ifndef M7_DATATYPE_H
#define M7_DATATYPE_H

#include <complex>

namespace datatype {
    /**
     * generic case: non-complex
     */
    template<typename T>
    struct is_complex_base_t : public std::false_type {};

    /**
     * complex case: is complex
     */
    template<typename T>
    struct is_complex_base_t<std::complex<T>> : public std::true_type {};

    template<typename T>
    struct remove_const_ref_t {
        typedef typename std::remove_reference<typename std::remove_const<T>::type>::type type;
    };

    /**
     * @return
     *  compile time constant bool: true if the template arg is a complex type
     */
    template<typename T>
    constexpr bool is_complex() {
        return is_complex_base_t<typename remove_const_ref_t<T>::type>::value;
    }

    /**
    * the std definition of "arithmetic" type only applies to primitive types, so here the definition is broadened to
    * include the complex numbers.
    */
    template<typename T>
    constexpr bool is_arithmetic() { return is_complex<T>() || std::is_arithmetic<T>::value; }


    template<typename T>
    static std::string name() { return typeid(T).name(); }

    template<>
    std::string name<char>() { return "char"; }

    template<>
    std::string name<short int>() { return "short int"; }

    template<>
    std::string name<int>() { return "int"; }

    template<>
    std::string name<long int>() { return "long int"; }

    template<>
    std::string name<long long int>() { return "long long int"; }

    template<>
    std::string name<unsigned char>() { return "unsigned char"; }

    template<>
    std::string name<unsigned short int>() { return "unsigned short int"; }

    template<>
    std::string name<unsigned int>() { return "unsigned short int"; }

    template<>
    std::string name<unsigned long int>() { return "unsigned long int"; }

    template<>
    std::string name<unsigned long long int>() { return "unsigned long long int"; }

    template<>
    std::string name<float>() { return "float"; }

    template<>
    std::string name<double>() { return "double"; }

    template<>
    std::string name<long double>() { return "long double"; }

    template<>
    std::string name<std::complex<float>>() { return "complex float"; }

    template<>
    std::string name<std::complex<double>>() { return "complex double"; }

    template<>
    std::string name<std::complex<long double>>() { return "complex long double"; }

    template<>
    std::string name<bool>() { return "bool"; }
}



#endif //M7_DATATYPE_H
