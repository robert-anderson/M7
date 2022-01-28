//
// Created by anderson on 1/27/22.
//

#ifndef M7_BUFFEREDFIELDS_H
#define M7_BUFFEREDFIELDS_H

#include "src/core/field/Fields.h"
#include "src/core/field/CompositeField.h"


struct BufferedFieldRow {
    Row m_internal_row;
    BufferedFieldRow():m_internal_row(){}
    BufferedFieldRow(const Row& row): m_internal_row(row){}
};

/**
 * Multiple inheritance is only used so that the m_internal_row member ctor is called before that of can be initialized
 * before that of the template arg
 */
template<typename T>
struct BufferedField : BufferedFieldRow, T {
    using T::operator=;
    Buffer m_internal_buffer;
    TableBase m_internal_table;

private:
    /**
     * contains common code for both ctors
     */
    void init() {
        m_internal_table.set_buffer(&m_internal_buffer);
        m_internal_row.m_table = &m_internal_table;
        m_internal_table.push_back();
        m_internal_row.restart();
    }

public:
    template<typename ...Args>
    BufferedField(Args&&... args):
        T(&m_internal_row, std::forward<Args>(args)...),
        m_internal_buffer("", 1), m_internal_table(m_internal_row.m_size) {
        init();
    }
//    BufferedField(const BufferedField& other):
//        BufferedFieldRow(other.m_internal_row), T(static_cast<const T&>(other)),
//        m_internal_buffer("", 1), m_internal_table(m_internal_row.m_size){
//        init();
//        // data is also copied:
//        *this = other;
//    }
//
//    BufferedField& operator=(const T& other) {
//        static_cast<T&>(*this) = other;
//        return *this;
//    }
};

namespace buffered {

    template<typename T, size_t nind>
    struct NdBitset : BufferedField<field::NdBitset<T, nind>> {
        typedef BufferedField<field::NdBitset<T, nind>> base_t;
        typedef typename field::NdBitset<T, nind>::inds_t inds_t;
        using field::NdBitset<T, nind>::operator=;
        explicit NdBitset(inds_t shape) : BufferedField<field::NdBitset<T, nind>>(shape){}
//        NdBitset(const field::NdBitset<T, nind>& field): BufferedField<field::NdBitset<T, nind>>(field){}
    };


    template<typename T>
    struct Bitset : NdBitset<T, 1ul> {
        explicit Bitset(size_t nbit): NdBitset<T, 1ul>({nbit}){}
    };


    template<typename T, size_t nind>
    struct Numbers : BufferedField<field::Numbers<T, nind>> {
        typedef BufferedField<field::Numbers<T, nind>> base_t;
        typedef typename field::Numbers<T, nind>::inds_t inds_t;
        using field::Numbers<T, nind>::operator=;
        explicit Numbers(inds_t shape) : base_t(shape){}
        Numbers(inds_t shape, T init_value) : base_t(shape){
            *this = init_value;
        }
//        explicit Numbers(const field::Numbers<T, nind>& field): BufferedField<field::Numbers<T, nind>>(field){}
//        Numbers(const Numbers& field): BufferedField<field::Numbers<T, nind>>(field){}
        Numbers& operator=(const Numbers& field){
            field::Numbers<T, nind>::operator=(field);
            return *this;
        }
    };

    template<typename T>
    struct Number : BufferedField<field::Number<T>> {
        using BufferedField<field::Number<T>>::operator=;
        Number(): BufferedField<field::Number<T>>({}){}
//        operator T&(){return (*this)[0];}
//        operator const T&() const {return (*this)[0];}
//        Number& operator=(const T& v){
//            static_cast<T&>(*this) = v;
//            return *this;
//        }
    };

    struct FrmOnv : BufferedField<field::FrmOnv> {
        using field::FrmOnv::operator=;
        FrmOnv(BasisDims bd) : BufferedField<field::FrmOnv>(bd){}
        FrmOnv(size_t nsite) : BufferedField<field::FrmOnv>(nsite){}
    };

    struct BosOnv : BufferedField<field::BosOnv> {
        using field::BosOnv::operator=;
        BosOnv(size_t nmode) : BufferedField<field::BosOnv>(nmode){}
        BosOnv(BasisDims bd) : BufferedField<field::BosOnv>(bd){}
    };

    struct FrmBosOnv : BufferedField<field::FrmBosOnv> {
        using field::FrmBosOnv::operator=;
        FrmBosOnv(BasisDims bd): BufferedField<field::FrmBosOnv>(bd){}
        FrmBosOnv(const field::FrmBosOnv& other): FrmBosOnv({other.m_frm.m_nsite, other.m_bos.m_nmode}){
            *this = other;
        }
        FrmBosOnv(const FrmBosOnv& other): FrmBosOnv(static_cast<const field::FrmBosOnv&>(other)){}
        FrmBosOnv& operator=(const FrmBosOnv& other){
            field::FrmBosOnv::operator=(other);
            return *this;
        }
    };

    typedef std::tuple<FrmOnv, FrmBosOnv, BosOnv> mbf_tup_t;

    template<size_t mbf_ind>
    using mbf_t = typename std::tuple_element<mbf_ind, mbf_tup_t>::type;
    typedef mbf_t<defs::mbf_type_ind> Mbf;

    struct MaeInds : BufferedField<field::MaeInds> {
        using field::MaeInds::operator=;
        using field::MaeInds::m_frm;
        using field::MaeInds::m_bos;
        MaeInds(size_t exsig): BufferedField<field::MaeInds>(exsig){}
    };

    struct SpecMomInds : BufferedField<field::SpecMomInds> {
        using field::SpecMomInds::operator=;
        SpecMomInds(size_t exsig): BufferedField<field::SpecMomInds>(exsig){}
    };
}

#endif //M7_BUFFEREDFIELDS_H
