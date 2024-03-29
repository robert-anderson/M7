//
// Created by Robert J. Anderson on 1/27/22.
//

#ifndef M7_BUFFEREDFIELDS_H
#define M7_BUFFEREDFIELDS_H

#include <M7_lib/field/Fields.h>
#include <M7_lib/field/CompositeField.h>

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

    template<typename ...Args>
    BufferedField(Args&&... args):
        T(&m_internal_row, std::forward<Args>(args)...), m_internal_table(m_internal_row.m_size) {
        m_internal_table.set_buffer(&m_internal_buffer);
        m_internal_row.m_table = &m_internal_table;
        m_internal_table.push_back();
        m_internal_row.restart();
    }
};

/**
 * the above-defined templated class is not intended to be directly usable: it only provides the data storage and a
 * generic forwarding constructor. The classes in this namespace inherit from BufferedField and re-implement the exact
 * ctors of the underlying Field or CompositeField along with the copy ctors and copy assignment operators.
 */
namespace buffered {

    template<typename T, uint_t nind>
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
        using BufferedField<field::Number<T>>::operator const T &;
        Number(): BufferedField<field::Number<T>>(){}
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

    template<typename T, uint_t nind>
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
        explicit Bitset(uint_t nbit): NdBitset<T, 1ul>({nbit}){}
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
        explicit FrmOnv(const sys::frm::Basis& basis) : BufferedField<field::FrmOnv>(basis){}
        FrmOnv(uint_t nsite) : FrmOnv(sys::frm::Basis(nsite)){}
        explicit FrmOnv(const sys::Basis& basis) : FrmOnv(basis.m_frm){}
        explicit FrmOnv(const sys::frm::Sector& sector) : FrmOnv(sector.m_basis){}
        explicit FrmOnv(const sys::Sector& sector) : FrmOnv(sector.m_frm){}

        FrmOnv(const FrmOnv& other) : FrmOnv(other.m_basis){
            *this = other;
        }
        FrmOnv(const field::FrmOnv& other) : FrmOnv(other.m_basis){
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
        explicit BosOnv(const sys::bos::Basis& basis) : BufferedField<field::BosOnv>(basis){}
        BosOnv(uint_t nmode, uint_t occ_cutoff=sys::bos::c_max_occ) : BosOnv(sys::bos::Basis(nmode, occ_cutoff)){}
        explicit BosOnv(const sys::Basis& basis) : BosOnv(basis.m_bos){}
        explicit BosOnv(const sys::bos::Sector& sector) : BosOnv(sector.m_basis){}
        explicit BosOnv(const sys::Sector& sector) : BosOnv(sector.m_bos){}

        BosOnv(const BosOnv& other) : BosOnv(other.m_basis){
            *this = other;
        }
        BosOnv(const field::BosOnv& other) : BosOnv(other.m_basis){
            *this = other;
        }
        BosOnv& operator=(const BosOnv& field){
            field::BosOnv::operator=(field);
            return *this;
        }
    };

    struct FrmBosOnv : BufferedField<field::FrmBosOnv> {
        using field::FrmBosOnv::operator=;
        explicit FrmBosOnv(const sys::frm::Basis& frm_basis, const sys::bos::Basis& bos_basis):
            BufferedField<field::FrmBosOnv>(frm_basis, bos_basis){}
        FrmBosOnv(uint_t nsite, uint_t nmode, uint_t bos_occ_cutoff=sys::bos::c_max_occ):
            FrmBosOnv(sys::frm::Basis(nsite), sys::bos::Basis(nmode, bos_occ_cutoff)){}
        explicit FrmBosOnv(const sys::Basis& basis): FrmBosOnv(basis.m_frm, basis.m_bos){}
        explicit FrmBosOnv(const sys::Sector& sector): FrmBosOnv(sector.basis()){}

        FrmBosOnv(const field::FrmBosOnv& other): FrmBosOnv(other.m_frm.m_basis, other.m_bos.m_basis){
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
        FrmXonv(const sys::frm::Basis& basis): BufferedField<field::FrmXonv>(basis){}
        FrmXonv(const sys::Basis& basis): BufferedField<field::FrmXonv>(basis){}
        FrmXonv(const sys::Sector& sector): BufferedField<field::FrmXonv>(sector){}

        FrmXonv(const field::FrmXonv& other): FrmXonv(other.m_ket.m_basis){
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
        BosXonv(const sys::bos::Basis& basis): BufferedField<field::BosXonv>(basis){}
        BosXonv(const sys::Basis& basis): BufferedField<field::BosXonv>(basis){}
        BosXonv(const sys::Sector& sector): BufferedField<field::BosXonv>(sector){}

        BosXonv(const field::BosXonv& other): BosXonv(other.m_ket.m_basis){
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
        FrmBosXonv(const sys::frm::Basis& frm_basis, const sys::bos::Basis& bos_basis):
                BufferedField<field::FrmBosXonv>(frm_basis, bos_basis){}
        FrmBosXonv(const sys::Basis& basis): BufferedField<field::FrmBosXonv>(basis){}
        FrmBosXonv(const sys::Sector& sector): BufferedField<field::FrmBosXonv>(sector){}

        FrmBosXonv(const field::FrmBosXonv& other): FrmBosXonv(other.m_ket.m_frm.m_basis, other.m_ket.m_bos.m_basis){
            *this = other;
        }
        FrmBosXonv(const FrmBosXonv& other): FrmBosXonv(static_cast<const field::FrmBosXonv&>(other)){}
        FrmBosXonv& operator=(const FrmBosXonv& other){
            field::FrmBosXonv::operator=(other);
            return *this;
        }
    };

    typedef std::tuple<FrmOnv, FrmBosOnv, BosOnv> mbf_tup_t;

    template<uint_t mbf_ind>
    using mbf_t = typename std::tuple_element<mbf_ind, mbf_tup_t>::type;
    typedef mbf_t<c_mbf_type_ind> Mbf;

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


    struct RdmInds : BufferedField<field::RdmInds> {
        using field::RdmInds::operator=;
        using field::RdmInds::m_frm;
        using field::RdmInds::m_bos;
        RdmInds(OpSig exsig): BufferedField<field::RdmInds>(exsig){}
        RdmInds(const RdmInds& other): BufferedField<field::RdmInds>(other.m_exsig){
            static_cast<field::RdmInds&>(*this) = other;
        }
    };

    struct SpecMomInds : BufferedField<field::SpecMomInds> {
        using field::SpecMomInds::operator=;
    };

    struct String : BufferedField<field::String> {
        using field::String::operator=;
        String(uint_t length): BufferedField<field::String>(length){}
    };
}

#endif //M7_BUFFEREDFIELDS_H
