//
// Created by rja on 23/06/22.
//

#ifndef M7_DATATYPE_H
#define M7_DATATYPE_H

#include <complex>

namespace dtype {
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

    /**
     * this is a recursive parameter pack expansion which casts every element to void, this is interpreted by the
     * compiler as the parameter being "used", and therefore suppresses the unused-parameter warnings when called, but
     * does not entail any excess overhead and should be optimized out
     */
    static void unused(){}

    template<typename T, typename ...Args>
    void unused(T&& first, Args&&... args) {
        (void)first;
        unused(args...);
    }

    template<typename T>
    constexpr T null() {
        return std::numeric_limits<T>::max();
    }
    template<typename T>
    constexpr T null(const T&) {
        return null<T>();
    }

    template<typename T>
    bool is_null(const T& v){
        return v==null(v);
    }

    template<typename T>
    static bool all_null(const T& v) {
        return is_null(v);
    }
    template<typename T, typename ...Args>
    bool all_null(const T& first, const Args&... rest){
        return is_null(first) && all_null(rest...);
    }

    template<typename T>
    static bool any_null(const T& v) {
        return is_null(v);
    }
    template<typename T, typename ...Args>
    bool any_null(const T& first, const Args&... rest){
        return is_null(first) || any_null(rest...);
    }

    template<typename T>
    void nullify(T& v){
        v = null(v);
    }

    template<typename T>
    void nullify(T* v){
        v = nullptr;
    }

    template<typename T, typename ...Args>
    void nullify(T& first, Args&... rest){
        nullify(first);
        nullify(rest...);
    }


    template<typename T, size_t n>
    static void nullify(std::array<T, n>& a) {
        for (uint i=0ul; i<n; ++i) nullify(a[i]);
    }
}



#endif //M7_DATATYPE_H
