//
// Created by rja on 09/02/2021.
//

#ifndef M7_FIELDSZ_H
#define M7_FIELDSZ_H

#include "NumberFieldZ.h"
#include "BitsetFieldZ.h"
#include "NdMultiFieldZ.h"

namespace fieldsz {

    template<typename T, size_t nind_item, size_t nind_element>
    struct NumberArrays : NdFieldZ<nind_item, NumberFieldZ<T, nind_item, nind_element>> {
        typedef NumberFieldZ<T, nind_item, nind_element> field_t;
        typedef NdFieldZ<nind_item, field_t> base_t;
        NumberArrays(RowZ* row, std::array<size_t, nind_item> item_shape, std::array<size_t, nind_element> element_shape) :
        base_t(row, item_shape, {element_shape}){}
    };

    template<typename T, size_t nind_element>
    struct NumberArray : NdFieldZ<0ul, NumberFieldZ<T, 0ul, nind_element>> {
        typedef NumberFieldZ<T, 0ul, nind_element> field_t;
        typedef NdFieldZ<0ul, field_t> base_t;
        NumberArray(RowZ* row, std::array<size_t, nind_element> element_shape) :
                base_t(row, {}, {element_shape}){}
    };

    template<typename T, size_t nind_item>
    struct Numbers : NdFieldZ<nind_item, NumberFieldZ<T, nind_item, 0ul>> {
        typedef NumberFieldZ<T, nind_item, 0ul> field_t;
        typedef NdFieldZ<nind_item, field_t> base_t;
        Numbers(RowZ* row, std::array<size_t, nind_item> item_shape) :
                base_t(row, item_shape, {}){}
    };

    template<typename T>
    struct Number : NdFieldZ<0ul, NumberFieldZ<T, 0ul, 0ul>> {
        typedef NumberFieldZ<T, 0ul, 0ul> field_t;
        typedef NdFieldZ<0ul, field_t> base_t;
        Number(RowZ* row) : base_t(row, {}, {}){}
    };

    template<size_t nind_item>
    struct FermionOnvs : NdFieldZ<nind_item, FermionOnvFieldZ<nind_item>> {
        typedef FermionOnvFieldZ<nind_item> field_t;
        typedef NdFieldZ<nind_item, field_t> base_t;
        FermionOnvs(RowZ* row, std::array<size_t, nind_item> item_shape, size_t nsite):
            base_t(row, item_shape, {nsite}){}
    };

    struct FermionOnv : NdFieldZ<0ul, FermionOnvFieldZ<0ul>> {
        typedef FermionOnvFieldZ<0ul> field_t;
        typedef NdFieldZ<0ul, field_t> base_t;
        FermionOnv(RowZ* row, size_t nsite):
                base_t(row, {}, {nsite}){}
    };

    template<size_t nind_item>
    using BosonOnvs = NdFieldZ<nind_item, BosonOnvFieldZ<nind_item>>;
    using BosonOnv = BosonOnvs<0ul>;


    template<size_t nind_item>
    struct FermiBosOnvs : NdMultiFieldZ<nind_item, FermionOnvFieldZ<nind_item>, BosonOnvFieldZ<nind_item>> {
        typedef NdMultiFieldZ<nind_item, FermionOnvFieldZ<nind_item>, BosonOnvFieldZ<nind_item>> base_t;
        FermiBosOnvs(RowZ *row, std::array<size_t, nind_item> item_shape, size_t nsite) :
                base_t(row, item_shape, {nsite}, {nsite}) {}

        const FermionOnvFieldZ<nind_item> &get_fonv() const {
            return this->template get<0>();
        }

        const BosonOnvFieldZ<nind_item> &get_bonv() const {
            return this->template get<1>();
        }

        FermionOnvFieldZ<nind_item> &get_fonv() {
            return this->template get<0>();
        }

        BosonOnvFieldZ<nind_item> &get_bonv() {
            return this->template get<1>();
        }
    };


    struct FermiBosOnv : NdMultiFieldZ<0ul, FermionOnvFieldZ<0ul>, BosonOnvFieldZ<0ul>> {
        typedef NdMultiFieldZ<0ul, FermionOnvFieldZ<0ul>, BosonOnvFieldZ<0ul>> base_t;
        FermiBosOnv(RowZ *row, size_t nsite) : base_t(row, {}, {nsite}, {{nsite}}) {}

        const FermionOnvFieldZ<0ul> &get_fonv() const {
            return this->template get<0>();
        }

        const BosonOnvFieldZ<0ul> &get_bonv() const {
            return this->template get<1>();
        }

        FermionOnvFieldZ<0ul> &get_fonv() {
            return this->template get<0>();
        }

        BosonOnvFieldZ<0ul> &get_bonv() {
            return this->template get<1>();
        }
    };


    template<size_t nind_item, bool enable_bosons=defs::enable_bosons>
    using Onvs = typename std::conditional<enable_bosons, FermiBosOnvs<nind_item>, FermionOnvs<nind_item>>::type;

    template<bool enable_bosons=defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;


    template<size_t nind_element>
    struct Flags : NdFieldZ<0ul, FlagFieldZ<void, nind_element>> {
        typedef FlagFieldZ<void, nind_element> field_t;
        typedef NdFieldZ<0ul, field_t> base_t;
        Flags(RowZ* row, std::array<size_t, nind_element> element_shape) : base_t(row, {}, {element_shape}){}
    };

    struct Flag : NdFieldZ<0ul, FlagFieldZ<void, 0ul>> {
        typedef FlagFieldZ<void, 0ul> field_t;
        typedef NdFieldZ<0ul, field_t> base_t;
        Flag(RowZ* row) : base_t(row, {}, {}){}
    };


#if 0
    template<typename T, size_t nind, size_t nind_element>
    struct NumberArrays : NdFieldZ<nind, NumberFieldZ<T, nind_element>> {
        NumberArrays(RowZ *row, std::array<size_t, nind> shape, std::array<size_t, nind_element> element_shape) :
                NdFieldZ<nind, NumberFieldZ<T, nind_element>>(row, shape, element_shape) {}
    };

    template<typename T, size_t nind_element>
    struct NumberArray : NdFieldZ<0ul, NumberFieldZ<T, nind_element>> {
        NumberArray(RowZ *row, std::array<size_t, nind_element> element_shape) :
                NdFieldZ<0ul, NumberFieldZ<T, nind_element>>(row, {}, element_shape) {}
    };

    template<typename T, size_t nind>
    struct Numbers : NdFieldZ<nind, NumberFieldZ<T, 0ul>> {
        Numbers(RowZ *row, std::array<size_t, nind> shape) :
                NdFieldZ<nind, NumberFieldZ<T, 0ul>>(row, shape, {}) {}

        template<typename ...Inds>
        const Numbers& select(Inds ...inds) const {
            NdFieldZ<nind, NumberFieldZ<T, 0ul>>::select_base(inds...);
            return *this;
        }

        operator const T &() const {
            return this->operator()();
        }

        operator T &() {
            return this->operator()();
        }
    };

    template<typename T>
    struct Number : Numbers<T, 0ul> {
        Number(RowZ *row) : Numbers<T, 0ul>(row, {}) {}
    };

    template<size_t nind>
    struct BosonOnvs : NdFieldZ<nind, BosonOnvFieldZ> {
        BosonOnvs(RowZ *row, std::array<size_t, nind> shape, size_t nmode) :
                NdFieldZ<nind, BosonOnvFieldZ>(row, shape, {nmode}) {}
    };

    struct BosonOnv : BosonOnvs<0ul> {
        BosonOnv(RowZ *row, size_t nmode) : BosonOnvs<0ul>(row, {}, nmode) {}
    };

    template<size_t nind_element>
    struct Flags : NdFieldZ<0ul, FlagFieldZ<nind_element>> {
        Flags(RowZ *row, std::array<size_t, nind_element> element_shape) :
                NdFieldZ<0ul, FlagFieldZ<nind_element>>(row, {}, element_shape) {}
    };

    struct Flag : Flags<0ul> {
        Flag(RowZ *row) : Flags<0ul>(row, {}) {}
    };

    template<size_t nind>
    struct FermionOnvs : NdFieldZ<nind, FermionOnvFieldZ> {
        FermionOnvs(RowZ *row, std::array<size_t, nind> shape, size_t nsite) :
                NdFieldZ<nind, FermionOnvFieldZ>(row, shape, {nsite}) {}
    };

    struct FermionOnv : FermionOnvs<0ul> {
        FermionOnv(RowZ *row, size_t nsite) : FermionOnvs<0ul>(row, {}, nsite) {}
    };

    template<size_t nind>
    struct FermiBosOnvs : NdMultiFieldZ<nind, FermionOnvFieldZ, BosonOnvFieldZ> {
        FermiBosOnvs(RowZ *row, std::array<size_t, nind> shape, size_t nsite) :
                NdMultiFieldZ<nind, FermionOnvFieldZ, BosonOnvFieldZ>(row, shape, {nsite}, {nsite}) {}

        using NdMultiFieldZ<nind, FermionOnvFieldZ, BosonOnvFieldZ>::get;

        const FermionOnvFieldZ &get_fonv() const {
            return this->template get<0>();
        }

        const FermionOnvFieldZ &get_bonv() const {
            return this->template get<1>();
        }
    };

    struct FermiBosOnv : FermiBosOnvs<0ul> {
        FermiBosOnv(RowZ *row, size_t nsite) : FermiBosOnvs<0ul>(row, {}, nsite) {}
    };

    template<size_t nind, bool enable_bosons = defs::enable_bosons>
    using Onvs = typename std::conditional<enable_bosons, FermiBosOnvs<nind>, FermionOnvs<nind>>::type;

    template<bool enable_bosons = defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;

#endif
}


#endif //M7_FIELDSZ_H
