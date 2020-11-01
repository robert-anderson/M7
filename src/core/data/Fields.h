//
// Created by RJA on 28/10/2020.
//

#ifndef M7_FIELDS_H
#define M7_FIELDS_H

#include "NdCompositeFieldSingle.h"
#include "NumericField.h"
#include "NumericArrayField.h"
#include "BosonOnvField.h"
#include "BitsetField.h"
#include "DeterminantField.h"
#include "ConfigurationField.h"

/*
 * some convenient wrappers for the field types
 */

namespace fields {

    /*
     * single composite fields...
     */

    template<typename T>
    struct Number : NdCompositeFieldSingle<NumericFieldX<T>, 0ul> {
        Number(TableX *table, std::string description) :
                NdCompositeFieldSingle<NumericFieldX<T>, 0ul>(table, {}, description, {}) {}
    };

    template<typename T, size_t nind>
    struct Numbers : NdCompositeFieldSingle<NumericFieldX<T>, nind> {
        template<typename...Args>
        Numbers(TableX *table, std::string description, Args...shape):
                NdCompositeFieldSingle<NumericFieldX<T>, nind>(table, {}, description, NdFormat<nind>(shape...)) {}
    };

    template<typename T, size_t nind_view>
    using NumberArray = NdCompositeFieldSingle<NumericArrayField<T, nind_view>, 0ul>;

    template<typename T, size_t nind_view, size_t nind_field>
    using NumberArrays = NdCompositeFieldSingle<NumericArrayField<T, nind_view>, nind_field>;

    struct Bitset : NdCompositeFieldSingle<BitsetFieldX, 0ul> {
        Bitset(TableX *table, size_t nbit, std::string description) :
                NdCompositeFieldSingle<BitsetFieldX, 0ul>(table, {nbit}, description, {}) {}
    };

    template<size_t nind>
    using Bitsets = NdCompositeFieldSingle<BitsetFieldX, nind>;

    struct Determinant : NdCompositeFieldSingle<DeterminantFieldX, 0ul> {
        Determinant(TableX *table, size_t nsite, std::string description) :
                NdCompositeFieldSingle<DeterminantFieldX, 0ul>(
                        table, {nsite}, description, {}) {}

        struct hash_fn {
            defs::hash_t operator()(const Determinant &composite, const_view_t v) const {
                return composite.m_field.hash(v);
            }
        };
    };

    template<size_t nind>
    using Determinants = NdCompositeFieldSingle<DeterminantFieldX, nind>;

    struct BosonOnv : NdCompositeFieldSingle<BosonOnvField, 0ul> {
        BosonOnv(TableX *table, size_t nmode, std::string description) :
                NdCompositeFieldSingle<BosonOnvField, 0ul>(
                        table, {nmode}, description, {}) {}
    };

    template<size_t nind>
    using BosonOnvs = NdCompositeFieldSingle<BosonOnvField, nind>;


    /*
     * true composite fields...
     */

    template<size_t nind, bool bosons>
    struct ConfigurationSelector {
    };

    template<size_t nind>
    struct ConfigurationSelector<nind, false> {
        typedef FermionField<nind> type;
    };

    template<size_t nind>
    struct ConfigurationSelector<nind, true> {
        typedef FermionBosonField<nind> type;
    };

    struct Configuration : ConfigurationSelector<0ul, defs::bosons>::type {
        Configuration(TableX *table, size_t nsite, size_t nmode, std::string description):
        ConfigurationSelector<0ul, defs::bosons>::type(table, nsite, nmode, description, {}){}

        struct hash_fn {
            defs::hash_t operator()(const Configuration &composite, const_view_t v) const {
                return composite.m_det.hash(v.m_det)^composite.m_perm.hash(v.m_perm);
            }
        };
    };

    template<size_t nind>
    using Configurations = typename ConfigurationSelector<nind, defs::bosons>::type;

}


#endif //M7_FIELDS_H
