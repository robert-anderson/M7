//
// Created by rja on 09/02/2021.
//

#ifndef M7_NUMBERFIELD_H
#define M7_NUMBERFIELD_H

#include "FieldBase.h"


struct NumberFieldBase : FieldBase {
    const size_t m_element_size, m_nelement;
    const bool m_is_complex;

    NumberFieldBase(Row* row, size_t element_size, size_t nelement, bool is_complex,
                    const std::type_info& type_info, std::string name=""):
    FieldBase(row, element_size * nelement, type_info, name),
    m_element_size(element_size), m_nelement(nelement), m_is_complex(is_complex){}

    virtual std::string stats_string() const = 0;

    virtual std::string format_string() const = 0;

    template<typename T>
    static std::string stats_string_element(const T& v){
        return utils::num_to_string(v);
    }

    template<typename T>
    static std::string stats_string_element(const std::complex<T>& v){
        return utils::num_to_string(v.real())+" "+utils::num_to_string(v.imag());
    }
};


namespace field_math {

    struct OperationBase {
    };

    struct Sqrt : OperationBase {
        template<typename T>
        T operator()(const T& v) const {
            return std::sqrt(v);
        }
    };

    struct Pow : OperationBase {
        const double m_exp;
        Pow(const double &exp) : m_exp(exp) {}
        template<typename T>
        T operator()(const T& v) const {
            return std::pow(v, m_exp);
        }
    };

    template<typename op_t, typename field_t>
    struct OpFieldPair {
        static_assert(std::is_base_of<OperationBase, op_t>::value, "Template arg must be derived from OperationBase");
        static_assert(std::is_base_of<NumberFieldBase, field_t>::value, "Template arg must be derived from NumberFieldBase");
        const field_t &m_rhs;
        const op_t m_op;

        template<typename ...Args>
        OpFieldPair(const field_t &rhs, Args &&... op_args): m_rhs(rhs), m_op(std::forward<Args>(op_args)...) {}

        template<typename T>
        T get(const size_t& ielement) const {
            return m_op(((const T*)(static_cast<const FieldBase&>(m_rhs).begin()))[ielement]);
        }
    };
}


template<typename T, size_t nind>
struct NdNumberField : NumberFieldBase {
    typedef const std::array<size_t, nind>& inds_t;
    const NdFormat<nind> m_format;

    const size_t& nelement() const {
        return m_format.nelement();
    }

    NdNumberField(Row* row, NdFormat<nind> format, std::string name=""):
            NumberFieldBase(row, sizeof(T), format.nelement(),
                            consts::is_complex<T>(), typeid(T), name), m_format(format){}

    NdNumberField(const NdNumberField& other):
            NdNumberField(other.row_of_copy(), other.m_format, other.m_name){}

    NdNumberField &operator=(const T &v) {
        std::fill((T*)begin(), (T*)end(), v);
        return *this;
    }

    NdNumberField &operator=(const std::vector<T> &v) {
        ASSERT(v.size()==nelement());
        std::memcpy(begin(), v.data(), m_size);
        return *this;
    }

    NdNumberField &operator=(const NdNumberField &v) {
        static_cast<FieldBase&>(*this) = v;
        return *this;
    }


    field_math::OpFieldPair<field_math::Sqrt, NdNumberField> sqrt() {
        return {*this};
    }

    template<typename field_t>
    field_math::OpFieldPair<field_math::Pow, field_t> pow(double exp) {
        return {*this, exp};
    }

    template<typename op_t>
    NdNumberField &operator=(const field_math::OpFieldPair<op_t, NdNumberField>& op_field_pair) {
        for (size_t ielement = 0ul; ielement<m_nelement; ++ielement)
            (*this)[ielement] = op_field_pair.template get<T>(ielement);
        return *this;
    }

    T sum() const {
        T tmp = 0;
        for (size_t ielement=0ul; ielement<m_nelement; ++ielement) tmp+=(*this)[ielement];
        return tmp;
    }

    T& operator[](const size_t& ielement) {
        return ((T *) begin())[ielement];
    }

    const T& operator[](const size_t& ielement) const {
        return ((const T *) begin())[ielement];
    }

    T& operator[](inds_t inds) {
        return ((T *) begin())[m_format.flatten(inds)];
    }

    const T& operator[](inds_t inds) const {
        return ((const T *) begin())[m_format.flatten(inds)];
    }

    std::string to_string() const override {
        std::string tmp;
        if (nind>0) tmp += "[";
        for (size_t ielement = 0ul; ielement<nelement(); ++ielement)
            tmp+=utils::num_to_string((*this)[ielement]) + " ";
        if (nind>0) tmp += "]";
        return tmp;
    }

    std::string stats_string() const override {
        std::string tmp;
        for (size_t ielement = 0ul; ielement<nelement(); ++ielement)
            tmp+= stats_string_element((*this)[ielement]);
        return tmp;
    }

    std::string format_string() const override {
        const std::string complex_dim_string = "real/imag (2)";
        if (!nind && m_is_complex) return complex_dim_string;
        return m_format.to_string() + (m_is_complex ? complex_dim_string : "");
    }

    defs::inds h5_shape() const override {
        return {m_format.shape().begin(), m_format.shape().end()};
    }

    std::vector<std::string> h5_dim_names() const override {
        if (!nind) return {};
        return {m_format.dim_names().begin(), m_format.dim_names().end()};
    }

    hid_t h5_type() const override {
        return hdf5::type<T>();
    }


};

template<typename T>
struct NumberField : NdNumberField<T, 0ul> {
    typedef NdNumberField<T, 0ul> base_t;
    using base_t::operator=;

    NumberField(Row* row, std::string name=""): base_t(row, {}, name){}

    operator T&() {
        return *(T*)FieldBase::begin();
    }

    operator const T&() const {
        return *(const T*)FieldBase::begin();
    }
};


#endif //M7_NUMBERFIELD_H
