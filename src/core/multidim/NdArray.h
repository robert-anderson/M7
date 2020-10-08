//
// Created by Robert John Anderson on 2020-02-22.
//

#ifndef M7_NDARRAY_H
#define M7_NDARRAY_H

#include <cstddef>

template <typename interface_t, size_t nind>
struct ArrayFormat {
    interface_t *m_interface;
    std::array<size_t, nind> m_shape;
    std::array<size_t, nind> m_strides;

    template<typename ...Args>
    ArrayFormat(interface_t *interface, Args&&... shape):m_interface(interface){
        static_assert(sizeof...(Args)>=nind, "not enough arguments to specify array shape");
        static_assert(sizeof...(Args)<=nind, "too many arguments to specify array shape");
        set_shape(std::forward<Args>(shape)...);
        m_strides.back() = 1ul;
        for (auto i = 2ul; i <= nind; i++) {
            m_strides[nind - i] = m_strides[nind - i + 1] * m_shape[nind - i + 1];
        }
    }

    void set_shape(){}
    template<typename ...Args>
    void set_shape(size_t first, Args... rest){
        m_shape[nind-sizeof...(rest)-1] = first;
        set_shape(rest...);
    }

    size_t partial_offset(size_t nind_unspec) const {return 0;}
    template<typename ...Args>
    size_t partial_offset(size_t nind_unspec, size_t first, Args... rest) const{
        return first*m_strides[nind-nind_unspec-1-sizeof...(Args)]+partial_offset(nind_unspec, rest...);
    }

    size_t nelement() const {
        return m_shape.front()*m_strides.front();
    }
};



template <typename interface_t, size_t nind, size_t nind_unspec>
struct ArrayIndexer {
    const ArrayFormat<interface_t, nind> &m_format;
    const size_t m_offset;
    static constexpr size_t nind_spec = nind-nind_unspec;

    ArrayIndexer(const ArrayFormat<interface_t, nind> &format) : m_format(format), m_offset(0ul) {}

    template<typename ...Args>
    ArrayIndexer(const ArrayIndexer<interface_t, nind, nind_unspec + sizeof...(Args)> &parent, Args... inds) :
            m_format(parent.m_format), m_offset(parent.m_offset+m_format.partial_offset(nind_unspec, inds...)){}

    size_t nelement() const {
        return m_format.m_shape[nind_spec]*m_format.m_strides[nind_spec];
    }

    template<typename T, size_t nind_unspec_after>
    struct Accessor{
        typedef ArrayIndexer<interface_t, nind, nind_unspec_after> type;
        template<typename ...Args>
        static type get(const ArrayIndexer<interface_t, nind, nind_unspec>& indexer, Args... inds){
            return ArrayIndexer<interface_t, nind, nind_unspec_after>(indexer, inds...);
        }
    };

    template<typename T>
    struct Accessor<T, 0ul>{
        typedef typename interface_t::View type;
        template<typename ...Args>
        static type get(const ArrayIndexer<interface_t, nind, nind_unspec>& indexer, Args... inds) {
            return indexer.m_format.m_interface->view(ArrayIndexer<interface_t, nind, 0ul>(indexer, inds...).m_offset);
        }
    };

    template<typename ...Args>
    typename Accessor<void, nind_unspec-sizeof...(Args)>::type
    operator()(Args... inds){
        return Accessor<void, nind_unspec-sizeof...(Args)>::get(*this, inds...);
    }
};

template<typename interface_t, size_t nind>
struct Array : public ArrayFormat<interface_t, nind> {
    template<typename ...Args>
    Array(interface_t *interface, Args&&... shape) :
            ArrayFormat<interface_t, nind>(interface, std::forward<Args>(shape)...){}

    ArrayIndexer<interface_t, nind, nind>
    operator()(){
        return ArrayIndexer<interface_t, nind, nind>(*this);
    }

    template<typename ...Args>
    typename ArrayIndexer<interface_t, nind, nind>::template Accessor<void, nind-sizeof...(Args)>::type
    operator()(Args... inds){
        return ArrayIndexer<interface_t, nind, nind>::template Accessor<void, nind-sizeof...(Args)>::get(*this, inds...);
    }
};

template<size_t nind>
struct TrivialArray : public Array<TrivialArray<nind>, nind> {
    template<typename ...Args>
    TrivialArray(Args&&... shape) :
            Array<TrivialArray<nind>, nind>(this, std::forward<Args>(shape)...){}

    struct View {
        size_t i;
        View(TrivialArray<nind>& array, size_t offset){
            i = offset;
        }
        operator size_t() const {
            return i;
        }
    };
    View view(const size_t& i){
        return View(*this, i);
    }
};


template<typename T, size_t nind>
struct NumericArray : public Array<NumericArray<T, nind>, nind> {

    using Array<NumericArray<T, nind>, nind>::nelement;

    std::vector<T> m_data;
    template<typename ...Args>
    NumericArray(Args&&... shape) :
            Array<NumericArray<T, nind>, nind>(this, std::forward<Args>(shape)...),
            m_data(nelement(), 0){}

    struct View {
        T& m_v;
        View(NumericArray<T, nind>& array, size_t offset):m_v(array.m_data[offset]){}

        operator T&() const {
            return m_v;
        }

        View& operator =(const T& v) {
            m_v = v;
        }
    };
    View view(const size_t& i){
        return View(*this, i);
    }
};


#endif //M7_NDARRAY_H
