//
// Created by rja on 28/06/22.
//

#ifndef M7_INTEGRALREADER_H
#define M7_INTEGRALREADER_H

#include "M7_lib/hdf5/File.h"
#include "M7_lib/io/FcidumpTextFileReader.h"
#include "M7_lib/linalg/Dense.h"

struct IntegralReader {
    struct IterData {
        uintv_t m_inds = uintv_t(4ul, ~0ul);
        ham_t m_value {};
        OpSig m_ranksig {};
        OpSig m_exsig {};

        /**
         * after setting m_inds, call this method to update the rank and excitation signatures
         */
        void update_ranksig_exsig();
    };

    virtual ~IntegralReader(){}
    /**
     * go to the next (inds, value) entry
     */
    virtual bool next(IterData& data) = 0;
    /**
     * go to the beginning of the 1e entries (needed for automatic permutational symmetry detection)
     */
    virtual void goto_first_1e() = 0;
    /**
     * go to the beginning of the 2e entries (needed for automatic permutational symmetry detection)
     */
    virtual void goto_first_2e() = 0;
    /**
     * get the core energy if the format exposes this value in a separate structure to the 1e, 2e entries
     */
    virtual ham_t ecore() const = 0;
    /**
     * whether the values data require a complex container
     */
    virtual bool complex_valued() const = 0;

    virtual bool spin_conserving(uint_t iex) const = 0;
};

struct CsvIntegralReader : IntegralReader {
    FcidumpTextFileReader m_reader;
    uint_t m_iline = 0ul;
    uint_t m_iline_first_1e = ~0ul;
    uint_t m_iline_first_2e = ~0ul;
    CsvIntegralReader(const FcidumpInfo& info);

    bool next(IterData& data) override;

    void goto_first_1e() override;

    void goto_first_2e() override;

    ham_t ecore() const override;

    bool complex_valued() const override;

    bool spin_conserving(uint_t iex) const override;
};

/**
 * assumes the existence of four file-level datasets:
 *  - 1 electron indices: (2, nentry) int64
 *  - 1 electron values: (nentry) double
 *  - 2 electron indices: (4, nentry) int64
 *  - 2 electron values: (nentry) double
 * and an attribute with the zero-electron contribution to the energy
 * only supports real integrals
 */
struct Hdf5IntegralReader : IntegralReader {
    struct KeyNames {
        const str_t m_ecore;
        const str_t m_1e_inds;
        const str_t m_1e_values;
        const str_t m_2e_inds;
        const str_t m_2e_values;
    };
private:
    hdf5::FileReader m_reader;
    const KeyNames m_names;
    dense::Matrix<int64_t> m_indices_2e;
    v_t<ham_comp_t> m_values_2e;
    dense::Matrix<int64_t> m_indices_1e;
    v_t<ham_comp_t> m_values_1e;
    uint_t m_iline = 0ul;
public:
    Hdf5IntegralReader(const FcidumpInfo& info, KeyNames names);
    bool next(IterData& data) override;

    void goto_first_1e() override;

    void goto_first_2e() override;

    ham_t ecore() const override;

    bool complex_valued() const override;

    bool spin_conserving(uint_t iex) const override;
};

/**
 * Molcas file-level keys are:
 * ['FOCK_INDEX', 'FOCK_VALUES', 'ORBITAL_ENERGIES', 'ORBITAL_INDEX', 'TWO_EL_INT_INDEX', 'TWO_EL_INT_VALUES']
 * note the unfortunate use of the misleading label "FOCK" in these keys, the quantity represented by that pair of
 * datasets is actually the 1-electron hamiltonian
 */
struct MolcasHdf5IntegralReader : Hdf5IntegralReader {
    MolcasHdf5IntegralReader(const FcidumpInfo& info);
};

#endif //M7_INTEGRALREADER_H
