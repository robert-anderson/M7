//
// Created by rja on 28/06/22.
//

#ifndef M7_INTEGRALREADER_H
#define M7_INTEGRALREADER_H

#include "M7_lib/hdf5/File.h"
#include "M7_lib/util/Exsig.h"
#include "M7_lib/io/FcidumpTextFileReader.h"
#include "M7_lib/linalg/Dense.h"

struct IntegralReader {
    struct IterData {
        uintv_t m_inds = uintv_t(4ul, ~0ul);
        ham_t m_value{};
        uint_t m_ranksig{};
        uint_t m_exsig{};

        void update_ranksig_exsig() {
            if (dtype::all_null(m_inds[0], m_inds[1], m_inds[2], m_inds[3])) m_ranksig = 0ul;
            else {
                REQUIRE_FALSE(dtype::any_null(m_inds[0], m_inds[1]),
                              "if not core energy entry, first two inds must be defined");
                m_ranksig = dtype::all_null(m_inds[2], m_inds[3]) ? exsig::ex_single : exsig::ex_double;
            }
            switch (m_ranksig) {
                case 0ul:
                    m_exsig = 0ul;
                    break;
                case exsig::ex_single:
                    m_exsig = m_inds[0] == m_inds[1] ? 0ul : exsig::ex_single;
                    break;
                case exsig::ex_double:
                    m_exsig = m_inds[0] == m_inds[2] ?
                              (m_inds[1] == m_inds[3] ? 0ul : exsig::ex_single) :
                              (m_inds[1] == m_inds[3] ? exsig::ex_single : exsig::ex_double);
                    break;
                default:
                    m_exsig = dtype::null(m_exsig);
            }
        }
    };

    virtual ~IntegralReader(){}
    virtual bool next(IterData& data) = 0;
    virtual void goto_first_1e() = 0;
    virtual void goto_first_2e() = 0;
    virtual ham_t ecore() const = 0;
    virtual bool complex_valued() const = 0;
};

struct CsvIntegralReader : IntegralReader {
    FcidumpTextFileReader m_reader;
    uint_t m_iline = 0ul;
    uint_t m_iline_first_1e = ~0ul;
    uint_t m_iline_first_2e = ~0ul;
    CsvIntegralReader(const FcidumpInfo& info, bool spin_major): m_reader(info.m_fname, spin_major){}
    bool next(IterData& data) override {
        auto out = m_reader.next(data.m_inds, data.m_value);
        if (!out) return false;
        data.m_ranksig = m_reader.ranksig(data.m_inds);
        data.m_exsig = m_reader.exsig(data.m_inds, data.m_ranksig);
        if (data.m_ranksig==exsig::ex_double && m_iline_first_2e==~0ul) m_iline_first_2e = m_iline;
        else if (data.m_ranksig==exsig::ex_single && m_iline_first_1e==~0ul) m_iline_first_1e = m_iline;
        ++m_iline;
        return out;
    }

    void goto_first_1e() override {
        REQUIRE_NE(m_iline_first_1e, ~0ul, "first 1-electron integral yet to be found");
        m_reader.reset(m_iline_first_1e);
        m_iline = m_iline_first_1e;
    }

    void goto_first_2e() override {
        REQUIRE_NE(m_iline_first_2e, ~0ul, "first 2-electron integral yet to be found");
        m_reader.reset(m_iline_first_2e);
        m_iline = m_iline_first_2e;
    }

    ham_t ecore() const override {
        return 0.0;
    }

    bool complex_valued() const override {
        return m_reader.m_complex_valued;
    }
};

/**
 * assumes the existence of four file-level datasets:
 *  - 1 electron indices: (2, nentry) int64
 *  - 1 electron values: (nentry) double
 *  - 2 electron indices: (4, nentry) int64
 *  - 2 electron values: (nentry) double
 * and an attribute with the zero-electron contribution to the energy
 */
struct Hdf5IntegralReader : IntegralReader {
    struct KeyNames {
        const std::string m_ecore;
        const std::string m_1e_inds;
        const std::string m_1e_values;
        const std::string m_2e_inds;
        const std::string m_2e_values;
    };
private:
    hdf5::FileReader m_reader;
    const KeyNames m_names;
    dense::Matrix<int64_t> m_indices_2e;
    std::vector<ham_t> m_values_2e;
    dense::Matrix<int64_t> m_indices_1e;
    std::vector<ham_t> m_values_1e;
    uint_t m_iline = 0ul;
public:
    Hdf5IntegralReader(const FcidumpInfo& info, KeyNames names, bool /*spin_major*/):
            m_reader(info.m_fname), m_names(std::move(names)),
            m_indices_2e(m_reader, m_names.m_2e_inds),
            m_values_2e(m_reader.read_data<std::vector<ham_t>>(m_names.m_2e_values)),
            m_indices_1e(m_reader, m_names.m_1e_inds),
            m_values_1e(m_reader.read_data<std::vector<ham_t>>(m_names.m_1e_values)){
        REQUIRE_EQ(m_indices_2e.nrow(), m_values_2e.size(),
                   "number of 2e matrix index arrays should match the number of values");
        REQUIRE_EQ(m_indices_1e.nrow(), m_values_1e.size(),
                   "number of 1e matrix index arrays should match the number of values");
    }
    bool next(IterData& data) override {
        if (m_iline== m_values_2e.size() + m_values_1e.size()) return false;
        if (m_iline < m_values_2e.size()) {
            const auto ientry = m_iline;
            m_indices_2e.get_row(ientry, data.m_inds);
            integer::shift(data.m_inds, false); // Fortran (1-based) to C (0-based) indexing
            data.update_ranksig_exsig();
            data.m_value = m_values_2e[ientry];
        }
        else {
            const auto ientry = m_iline-m_values_2e.size();
            data.m_inds[0] = m_indices_1e(ientry, 0) - 1;
            data.m_inds[1] = m_indices_1e(ientry, 1) - 1;
            dtype::nullify(data.m_inds[2], data.m_inds[3]);
            data.update_ranksig_exsig();
            data.m_value = m_values_1e[ientry];
        }
        ++m_iline;
        return true;
    }

    void goto_first_1e() override {
        m_iline = m_values_2e.size();
    }

    void goto_first_2e() override {
        m_iline = 0;
    }

    ham_t ecore() const override {
        return m_reader.read_attr<double>(m_names.m_ecore);
    }

    bool complex_valued() const override {
        return false;
    }
};

/**
 * Molcas file-level keys are:
 * ['FOCK_INDEX', 'FOCK_VALUES', 'ORBITAL_ENERGIES', 'ORBITAL_INDEX', 'TWO_EL_INT_INDEX', 'TWO_EL_INT_VALUES']
 * note the unfortunate use of the misleading label "FOCK" in these keys, the quantity represented by that pair of
 * datasets is actually the 1-electron hamiltonian
 */
struct MolcasHdf5IntegralReader : Hdf5IntegralReader {
    MolcasHdf5IntegralReader(const FcidumpInfo& info, bool spin_major):
        Hdf5IntegralReader(info,
               {"CORE_ENERGY",
                "FOCK_INDEX",
                "FOCK_VALUES",
                "TWO_EL_INT_INDEX",
                "TWO_EL_INT_VALUES"}, spin_major){}
};

#endif //M7_INTEGRALREADER_H
