//
// Created by anderson on 1/27/22.
//

#ifndef M7_BUFFEREDFIELDS2_H
#define M7_BUFFEREDFIELDS2_H

struct BufferedFieldRow {
    Row m_internal_row;
    BufferedFieldRow():m_internal_row(){}
    BufferedFieldRow(const Row& row): m_internal_row(row){}
};

/**
 * Multiple inheritance is only used so that the m_internal_row member ctor is called before that
 * of can be initialized before that of the template arg
 */
template<typename T>
struct BufferedField2 : BufferedFieldRow, T {
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
    BufferedField2(Args&&... args): T(&m_internal_row, std::forward<Args>(args)...),
        m_internal_buffer("", 1), m_internal_table(m_internal_row.m_size) {
        init();
    }

    BufferedField2(const BufferedField2& other): Wrapp(other.m_internal_row), T(static_cast<const T&>(other)),
        m_internal_buffer("", 1), m_internal_table(m_internal_row.m_size){
        init();
        // data is also copied:
        *this = other;
    }

    BufferedField2& operator=(const T& other) {
        static_cast<T&>(*this) = other;
        return *this;
    }

    BufferedField2& operator=(const BufferedField2& other) {
        static_cast<T&>(*this) = static_cast<const T&>(other);
        return *this;
    }
};

namespace buffered2 {

    template<typename T, size_t nind>
    struct NdBitset : BufferedField2<field::NdBitset<T, nind>> {
        typedef BufferedField2<field::NdBitset<T, nind>> base_t;
        typedef typename field::NdBitset<T, nind>::inds_t inds_t;
        using field::NdBitset<T, nind>::operator=;
        NdBitset(inds_t shape) : BufferedField2<field::NdBitset<T, nind>>({nullptr, shape}){}
        NdBitset(const field::NdBitset<T, nind>& field): BufferedField2<field::NdBitset<T, nind>>(field){}
    };


    template<typename T>
    struct Bitset : NdBitset<T, 1ul> {
        Bitset(size_t nbit): NdBitset<T, 1ul>({nbit}){}
    };


    template<typename T, size_t nind>
    struct Numbers : BufferedField<field::Numbers<T, nind>> {
        typedef BufferedField<field::Numbers<T, nind>> base_t;
        typedef typename field::Numbers<T, nind>::inds_t inds_t;
        using field::Numbers<T, nind>::operator=;
        Numbers(inds_t shape) : base_t({nullptr, shape}){}
        Numbers(inds_t shape, T init_value) : base_t({nullptr, shape}){
            *this = init_value;
        }
        Numbers(const field::Numbers<T, nind>& field): BufferedField<field::Numbers<T, nind>>(field){}
        Numbers(const Numbers& field): BufferedField<field::Numbers<T, nind>>(field){}
        Numbers& operator=(const Numbers& field){
            base_t::operator=(field);
            return *this;
        }
    };

    template<typename T>
    struct Number : BufferedField<field::Number<T>> {
        using BufferedField<field::Number<T>>::operator=;
        Number(): BufferedField<field::Number<T>>({{}}){}
        operator T&(){return (*this)[0];}
        operator const T&() const {return (*this)[0];}
        Number& operator=(const T& v){
            static_cast<T&>(*this) = v;
            return *this;
        }
    };

    struct FrmOnv : BufferedField<field::FrmOnv> {
        using field::FrmOnv::operator=;
        FrmOnv(BasisDims bd) : BufferedField<field::FrmOnv>({nullptr, bd}){}
        FrmOnv(size_t nsite) : BufferedField<field::FrmOnv>({nullptr, nsite}){}
        FrmOnv& operator=(const FrmOnv& other){
            base_t::operator=(other);
            return *this;
        }
        FrmOnv(const FrmOnv& other): FrmOnv({other.m_nsite, 0ul}){}
    };

    struct BosOnv : BufferedField<field::BosOnv> {
        using field::BosOnv::operator=;
        BosOnv(size_t nmode) : BufferedField<field::BosOnv>({nullptr, nmode}){}
        BosOnv(BasisDims bd) : BufferedField<field::BosOnv>({nullptr, bd}){}
    };


    struct FrmBosOnv : BufferedMultiField<field::FrmBosOnv> {
        using field::FrmBosOnv::operator=;
        FrmBosOnv(BasisDims bd):
                BufferedMultiField<field::FrmBosOnv>({nullptr, bd}){}
    };

    typedef std::tuple<FrmOnv, FrmBosOnv, BosOnv> mbf_tup_t;

    template<size_t mbf_ind>
    using mbf_t = typename std::tuple_element<mbf_ind, mbf_tup_t>::type;
    typedef mbf_t<defs::mbf_type_ind> Mbf;

    struct MaeInds : BufferedField<field::MaeInds> {
        using field::MaeInds::operator=;
        using field::MaeInds::m_frm;
        using field::MaeInds::m_bos;
        MaeInds(size_t exsig):
        BufferedField<field::MaeInds>({nullptr, exsig}){}
    };

    struct SpecMomInds : BufferedMultiField<field::SpecMomInds> {
        using field::SpecMomInds::operator=;
        SpecMomInds(size_t exsig):
        BufferedMultiField<field::SpecMomInds>({nullptr, exsig}){}
    };
}

template<typename field_t>
struct SingleFieldRow : Row {
    field_t m_field;
    template<typename ...Args>
    SingleFieldRow(Args... args): Row(), m_field(this, args...){}

    field_t &key_field() {
        return m_field;
    };

};

#endif //M7_BUFFEREDFIELDS2_H
