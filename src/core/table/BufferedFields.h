//
// Created by anderson on 1/27/22.
//

#ifndef M7_BUFFEREDFIELDS_H
#define M7_BUFFEREDFIELDS_H

#include "field/Fields.h"
#include "field/CompositeField.h"

struct BufferedFieldRow {
    Row m_internal_row;
    Buffer m_internal_buffer;
    BufferedFieldRow();
    BufferedFieldRow(const Row& row);
};

/**
 * BufferedFields allow for the use of Fields and CompositeFields without referencing an external Table to provide the
 * data buffer. Instead this class in effect provides subclasses of Fields and CompositeFields which own their own data
 * buffers.
 *
 * Multiple inheritance is only used so that the m_internal_row member ctor is called before that of can be initialized
 * before that of the template arg
 */
template<typename T>
struct BufferedField : BufferedFieldRow, T {
    using T::operator=;
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
        T(&m_internal_row, std::forward<Args>(args)...), m_internal_table(m_internal_row.m_size) {
        init();
    }
};

/**
 * the above-defined templated class is not intended to be directly usable: it only provides the data storage and a
 * generic forwarding constructor. The classes in this namespace inherit from BufferedField and re-implement the exact
 * ctors of the underlying Field or CompositeField along with the copy ctors and copy assignment operators.
 */
namespace buffered {

    template<typename T, size_t nind>
    struct Numbers : BufferedField<field::Numbers<T, nind>> {
        typedef BufferedField<field::Numbers<T, nind>> base_t;
        typedef typename field::Numbers<T, nind>::inds_t inds_t;
        using field::Numbers<T, nind>::operator=;
        explicit Numbers(inds_t shape) : base_t(shape){}
        Numbers(inds_t shape, T init_value) : base_t(shape){
            *this = init_value;
        }
        Numbers(const Numbers& other): Numbers(other.m_format.m_shape){
            *this = other;
        }
        Numbers(const field::Numbers<T, nind>& other): Numbers(other.m_format.m_shape){
            *this = other;
        }
        Numbers& operator=(const Numbers& field){
            field::Numbers<T, nind>::operator=(field);
            return *this;
        }
    };

    template<typename T>
    struct Number : BufferedField<field::Number<T>> {
        using BufferedField<field::Number<T>>::operator=;
        using BufferedField<field::Number<T>>::operator==;
        using BufferedField<field::Number<T>>::operator T &;
        Number(): BufferedField<field::Number<T>>({}){}
        Number(const Number& other): Number() {
            *this = other;
        }
        Number(const field::Number<T>& other): Number() {
            *this = other;
        }
        Number& operator=(const Number<T>& other) {
            field::Number<T>::operator=(other);
            return *this;
        }
    };

    template<typename T, size_t nind>
    struct NdBitset : BufferedField<field::NdBitset<T, nind>> {
        typedef BufferedField<field::NdBitset<T, nind>> base_t;
        typedef typename field::NdBitset<T, nind>::inds_t inds_t;
        using field::NdBitset<T, nind>::operator=;
        explicit NdBitset(inds_t shape) : BufferedField<field::NdBitset<T, nind>>(shape){}
        NdBitset(const NdBitset& other) : NdBitset(other.m_format.m_shape){
            *this = other;
        }
        NdBitset(const field::NdBitset<T, nind>& other) : NdBitset(other.m_format.m_shape){
            *this = other;
        }
        NdBitset& operator=(const NdBitset& field){
            field::NdBitset<T, nind>::operator=(field);
            return *this;
        }
    };

    template<typename T>
    struct Bitset : NdBitset<T, 1ul> {
        using field::Bitset<T>::operator=;
        explicit Bitset(size_t nbit): NdBitset<T, 1ul>({nbit}){}
        Bitset(const Bitset& other) : Bitset(other.m_format.m_nelement){
            *this = other;
        }
        Bitset(const field::Bitset<T>& other) : Bitset(other.m_format.m_nelement){
            *this = other;
        }
        Bitset& operator=(const Bitset& field){
            field::Bitset<T>::operator=(field);
            return *this;
        }
    };

    struct FrmOnv : BufferedField<field::FrmOnv> {
        using field::FrmOnv::operator=;
        FrmOnv(BasisData bd) : BufferedField<field::FrmOnv>(bd){}
        FrmOnv(size_t nsite) : BufferedField<field::FrmOnv>(nsite){}

        FrmOnv(const FrmOnv& other) : FrmOnv(other.m_nsite){
            *this = other;
        }
        FrmOnv(const field::FrmOnv& other) : FrmOnv(other.m_nsite){
            *this = other;
        }
        FrmOnv& operator=(const FrmOnv& field){
            field::FrmOnv::operator=(field);
            return *this;
        }
    };

    struct BosOnv : BufferedField<field::BosOnv> {
        using field::BosOnv::operator=;
        using field::BosOnv::operator==;
        BosOnv(size_t nmode) : BufferedField<field::BosOnv>(nmode){}
        BosOnv(BasisData bd) : BufferedField<field::BosOnv>(bd){}
        BosOnv(const BosOnv& other) : BosOnv(other.m_nmode){
            *this = other;
        }
        BosOnv(const field::BosOnv& other) : BosOnv(other.m_nmode){
            *this = other;
        }
        BosOnv& operator=(const BosOnv& field){
            field::BosOnv::operator=(field);
            return *this;
        }
    };

    struct FrmBosOnv : BufferedField<field::FrmBosOnv> {
        using field::FrmBosOnv::operator=;
        FrmBosOnv(BasisData bd): BufferedField<field::FrmBosOnv>(bd){}
        FrmBosOnv(const field::FrmBosOnv& other): FrmBosOnv({other.m_frm.m_nsite, other.m_bos.m_nmode}){
            *this = other;
        }
        FrmBosOnv(const FrmBosOnv& other): FrmBosOnv(static_cast<const field::FrmBosOnv&>(other)){}
        FrmBosOnv& operator=(const FrmBosOnv& other){
            field::FrmBosOnv::operator=(other);
            return *this;
        }
    };

    struct FrmXonv : BufferedField<field::FrmXonv> {
        using field::FrmXonv::operator=;
        FrmXonv(size_t nsite): BufferedField<field::FrmXonv>(nsite){}
        FrmXonv(BasisData bd): FrmXonv(bd.m_nsite){}
        FrmXonv(const field::FrmXonv& other): FrmXonv(other.m_ket.m_nsite){
            *this = other;
        }
        FrmXonv(const FrmXonv& other): FrmXonv(static_cast<const field::FrmXonv&>(other)){}
        FrmXonv& operator=(const FrmXonv& other){
            field::FrmXonv::operator=(other);
            return *this;
        }
    };

    struct BosXonv : BufferedField<field::BosXonv> {
        using field::BosXonv::operator=;
        BosXonv(size_t nmode): BufferedField<field::BosXonv>(nmode){}
        BosXonv(BasisData bd): BosXonv(bd.m_nmode){}
        BosXonv(const field::BosXonv& other): BosXonv(other.m_ket.m_nmode){
            *this = other;
        }
        BosXonv(const BosXonv& other): BosXonv(static_cast<const field::BosXonv&>(other)){}
        BosXonv& operator=(const BosXonv& other){
            field::BosXonv::operator=(other);
            return *this;
        }
    };


    struct FrmBosXonv : BufferedField<field::FrmBosXonv> {
        using field::FrmBosXonv::operator=;
        FrmBosXonv(BasisData bd): BufferedField<field::FrmBosXonv>(bd){}
        FrmBosXonv(const field::FrmBosXonv& other):
                FrmBosXonv({other.m_ket.m_frm.m_nsite, other.m_ket.m_bos.m_nmode}){
            *this = other;
        }
        FrmBosXonv(const FrmBosXonv& other): FrmBosXonv(static_cast<const field::FrmBosXonv&>(other)){}
        FrmBosXonv& operator=(const FrmBosXonv& other){
            field::FrmBosXonv::operator=(other);
            return *this;
        }
    };

    typedef std::tuple<FrmOnv, FrmBosOnv, BosOnv> mbf_tup_t;

    template<size_t mbf_ind>
    using mbf_t = typename std::tuple_element<mbf_ind, mbf_tup_t>::type;
    typedef mbf_t<defs::mbf_type_ind> Mbf;

    template<typename T=void>
    struct selector {
        typedef void type;
    };
    template<>
    struct selector<FrmOnv> {
        typedef FrmOnvField type;
    };
    template<>
    struct selector<FrmBosOnv> {
        typedef FrmBosOnvField type;
    };
    template<>
    struct selector<BosOnv> {
        typedef BosOnvField type;
    };

    template<typename T>
    using from_field_t = typename selector<T>::type;


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
