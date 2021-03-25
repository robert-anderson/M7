//
// Created by rja on 05/11/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <type_traits>
#include <src/defs.h>
#include "FermiBosHamiltonian.h"
#include "src/core/nd/NdArray.h"

template<bool enable_bosons = defs::enable_bosons>
using Hamiltonian = typename std::conditional<enable_bosons, FermiBosHamiltonian, FermionHamiltonian>::type;

#if 0

/*
 * every ...HamData type must define the following methods:
 *
 * get_coeff_0000
 * get_coeff_1100
 * get_coeff_2200
 * get_coeff_1101
 * get_coeff_1110
 * get_coeff_1111
 *
 * more can be added later, for example for higher-order interactions, but this set covers
 * abinitio, hubbard and hubbard-holstein.
 *
 * these are not implemented by inheritance to avoid runtime polymorphic dispatch overheads
 */


struct FermionIntegralArrays {
    const size_t m_nelec;
    const size_t m_nsite;
    const bool m_spin_conserving_1e, m_spin_conserving_2e;
    const bool m_complex_valued;
    const size_t m_int_2e_rank;
    defs::ham_t m_int_0;
    typedef Integrals_1e<defs::ham_t, defs::isym_1e> ints1_t;
    typedef Integrals_2e<defs::ham_t, defs::isym_2e> ints2_t;
    ints1_t m_int_1;
    ints2_t m_int_2;

    FermionIntegralArrays(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e,
                          bool spin_conserving_2e, bool complex_valued, bool spin_resolved, size_t int_2e_rank) :
            m_nelec(nelec), m_nsite(nsite),
            m_spin_conserving_1e(spin_conserving_1e),
            m_spin_conserving_2e(spin_conserving_2e),
            m_complex_valued(complex_valued),
            m_int_2e_rank(int_2e_rank),
            m_int_1(nsite, spin_resolved),
            m_int_2(nsite, spin_resolved) {}

    FermionIntegralArrays(const FcidumpFileReader &file_reader) :
            FermionIntegralArrays(file_reader.nelec(), file_reader.nspatorb(),
                                  file_reader.spin_conserving_1e(),
                                  file_reader.spin_conserving_2e(),
                                  file_reader.m_complex_valued,
                                  file_reader.spin_resolved(),
                                  file_reader.int_2e_rank()) {
        defs::inds inds(4);
        defs::ham_t value;

        log::info("Reading fermion Hamiltonian from FCIDUMP file \"" + file_reader.m_fname + "\"...");
        log::info("Loading fermion Hamiltonian from FCIDUMP...");
        while (file_reader.next(inds, value)) {
            if (ints2_t::valid_inds(inds)) m_int_2.set(inds, value);
            else if (ints1_t::valid_inds(inds)) m_int_1.set(inds, value);
            else if (inds[0] == ~0ul) m_int_0 = value;
        }
        mpi::barrier();
        log::info("FCIDUMP loading complete.");
    }

    FermionIntegralArrays(const std::string &fname, bool spin_major) :
            FermionIntegralArrays(FcidumpFileReader(fname, spin_major)) {}

};



struct FermionConnection {
    defs::inds m_ann;
    defs::inds m_cre;
    defs::inds m_com;
    bool m_phase;

    FermionConnection(size_t nsite) {
        m_ann.reserve(2 * nsite);
        m_cre.reserve(2 * nsite);
        m_com.reserve(2 * nsite);
    }

    size_t nann() const {
        return m_ann.size();
    }

    size_t ncre() const {
        return m_cre.size();
    }

    size_t ncom() const {
        return m_com.size();
    }

    void zero() {
        m_cre.clear();
        m_ann.clear();
        m_com.clear();
    }

    void add_cre(const size_t &i) {
        ASSERT(m_cre.size() < m_cre.capacity());
        m_cre.push_back(i);
    }

    void add_ann(const size_t &i) {
        ASSERT(m_ann.size() < m_ann.capacity());
        m_ann.push_back(i);
    }

    void add(const size_t &ann, const size_t &cre) {
        add_ann(ann);
        add_cre(cre);
    }

    void add(const size_t &ann1, const size_t &ann2, const size_t &cre1, const size_t &cre2) {
        add_ann(ann1);
        add_ann(ann2);
        add_cre(cre1);
        add_cre(cre2);
    }

    void sort() {
        std::sort(m_cre.begin(), m_cre.end());
        std::sort(m_ann.begin(), m_ann.end());
    }

    void connect(const fieldsx::FermionOnv &in, const fieldsx::FermionOnv &out) {
        ASSERT(!in.is_zero())
        zero();

        size_t nperm = 0ul;
        auto nann_found = 0ul;
        auto ncre_found = 0ul;

        defs::data_t in_work, out_work, work;
        for (size_t idataword = 0ul; idataword < in.m_item_dsize; ++idataword) {
            in_work = in.get_dataword(idataword);
            out_work = out.get_dataword(idataword);
            work = in_work & ~out_work;
            while (work) add_ann(bit_utils::next_setbit(work) + idataword * defs::nbit_data);
            work = out_work & ~in_work;
            while (work) add_cre(bit_utils::next_setbit(work) + idataword * defs::nbit_data);
            work = in_work & out_work;
            while (work) {
                auto setbit = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
                while (nann_found != nann() && m_ann[nann_found] < setbit) {
                    // an annihilation operator has been passed in the iteration over common indices
                    ++nann_found;
                    nperm += ncom();
                }
                while (ncre_found != ncre() && m_cre[ncre_found] < setbit) {
                    // a creation operator has been passed in the iteration over common indices
                    ++ncre_found;
                    nperm += ncom();
                }
                m_com.push_back(setbit);
            }
        }
        while (nann_found < nann()) {
            ++nann_found;
            nperm += ncom();
        }
        while (ncre_found < ncre()) {
            ++ncre_found;
            nperm += ncom();
        }
        m_phase = nperm & 1ul;
    }

    /*
     * compute phase of connection applied to "in" ONV
     */
    void apply(const fields::Onv<0> &in) {
        ASSERT(std::is_sorted(m_cre.begin(), m_cre.end()));
        ASSERT(std::is_sorted(m_ann.begin(), m_ann.end()));
        m_com.clear();
        size_t nperm = 0ul;

        auto ann_iter = m_ann.begin();
        auto cre_iter = m_cre.begin();

        for (size_t idataword = 0ul; idataword < in.m_item_dsize; ++idataword) {
            auto work = in.get_dataword(idataword);
            while (work) {
                auto setbit = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
                if (ann_iter != m_ann.end() && setbit == *ann_iter) {
                    ann_iter++;
                    nperm += ncom();
                    continue;
                }
                // check we aren't trying to create an electron in an occupied orbital
                ASSERT((cre_iter == m_cre.end()) || (setbit != *cre_iter));
                while (cre_iter != m_cre.end() && *cre_iter < setbit) {
                    cre_iter++;
                    nperm += ncom();
                }
                m_com.push_back(setbit);
            }
        }
        while (cre_iter != m_cre.end()) {
            cre_iter++;
            nperm += ncom();
        }
        m_phase = nperm & 1ul;
    }


    void apply(const fields::Onv<0> &in, fields::Onv<0> &out) {
        apply(in);
        ASSERT(!in.is_zero());
#ifndef NDEBUG
        for (const auto &i: m_ann) ASSERT(in.get(i));
        for (const auto &i: m_cre) ASSERT(!in.get(i));
#endif
        out = in;
        for (const auto &i: m_ann) out.clr(i);
        for (const auto &i: m_cre) out.set(i);
        ASSERT(in.nsetbit() == out.nsetbit());
    }
};

struct BosonConnection {
    defs::inds m_ann;
    defs::inds m_cre;
    defs::inds m_com;

    size_t nann() const {
        return m_ann.size();
    }

    size_t ncre() const {
        return m_cre.size();
    }

    size_t ncom() const {
        return m_com.size();
    }

};

struct ConnectionGeneral {
    FermionConnection m_fermion;
    BosonConnection m_boson;

    ConnectionGeneral(size_t nsite) : m_fermion(nsite) {

    }

    void connect(const fieldsx::Onv &in, const fieldsx::Onv &out) {
        m_fermion.connect(in.m_fonv, out.m_fonv);
    }

    void apply(const fieldsx::Onv &in, fieldsx::Onv &out) {
        m_fermion.apply(in.m_fonv, out.m_fonv);
    }

    bool phase() const {
        return m_fermion.m_phase;
    }

    static constexpr size_t c_max_nop_bits = 2;
    static constexpr size_t c_max_nop = 1ul << c_max_nop_bits;

    static constexpr size_t signature(size_t ncre, size_t nann) {
        return ncre | nann << c_max_nop_bits;
    }

    static constexpr size_t signature(size_t ncre_f, size_t nann_f, size_t ncre_b, size_t nann_b) {
        return signature(ncre_f, nann_f) | signature(ncre_b, nann_b) << (2 * c_max_nop_bits);
    }

    size_t get_signature() const {
        ASSERT(m_fermion.ncre() < c_max_nop);
        ASSERT(m_fermion.ncre() < c_max_nop);
        if (defs::enable_bosons) {
            ASSERT(m_boson.ncre() < c_max_nop);
            ASSERT(m_boson.ncre() < c_max_nop);
            return signature(m_fermion.ncre(), m_fermion.nann(), m_boson.ncre(), m_boson.nann());
        }
        return signature(m_fermion.ncre(), m_fermion.nann());
    }
};

#if 0
struct HamDataX {
};

struct HamData {
    /*
     * min and max excitation levels corresponding to each term:
     */
    static constexpr size_t c_nsignature = 1ul << (ConnectionGeneral::c_max_nop*4);
    typedef std::vector<std::pair<size_t, size_t>> insert_pairs_t;
    const std::array<insert_pairs_t, c_nsignature> m_promotions;

    std::array<insert_pairs_t, c_nsignature> make_promotions()

    static constexpr size_t c_max_nop = 2;
    const NdFormat<4> m_exlvl_format;
    const size_t m_nsite;
    typedef std::pair<std::array<size_t, 4>, std::array<size_t, 4>> exlvl_pair_t;
    const exlvl_pair_t m_exlvl_1100; // c+ c
    const exlvl_pair_t m_exlvl_2200; // c+ c+ c c
    const exlvl_pair_t m_exlvl_1101; // c+ c b
    const exlvl_pair_t m_exlvl_1110; // c+ c b+
    const exlvl_pair_t m_exlvl_0011; // b+ b

    HamData(size_t nsite, exlvl_pair_t exlvl_1100, exlvl_pair_t exlvl_2200, exlvl_pair_t exlvl_1101,
            exlvl_pair_t exlvl_1110, exlvl_pair_t exlvl_0011) : m_exlvl_format(c_max_nop),
                                                                m_nsite(nsite), m_exlvl_1100(exlvl_1100),
                                                                m_exlvl_2200(exlvl_2200), m_exlvl_1101(exlvl_1101),
                                                                m_exlvl_1110(exlvl_1110), m_exlvl_0011(exlvl_0011) {}

    static exlvl_pair_t conserve_fermion_exlvl(size_t min, size_t max) {
        ASSERT(min <= max);
        return {{min, min, 0, 0},
                {max, max, 0, 0}};
    }

    static exlvl_pair_t null_exlvl() {
        return {{0, 0, 0, 0},
                {0, 0, 0, 0}};
    }

};

struct AbinitioHamData : HamData {
    FermionIntegralArrays m_ints;

    AbinitioHamData(FermionIntegralArrays &&ints) :
            HamData(ints.m_nsite,
                    conserve_fermion_exlvl(0, 1),
                    conserve_fermion_exlvl(0, 2),
                    null_exlvl(), null_exlvl(), null_exlvl()), m_ints(std::move(ints)) {}

    defs::ham_t get_coeff_0000() const {
        return m_ints.m_int_0;
    }

    defs::ham_t get_coeff_1100(const size_t& i, const size_t& j) const {
        return m_ints.m_int_1.get(i, j);
    }

    defs::ham_t get_coeff_2200(const size_t& i, const size_t& j, const size_t& k, const size_t& l) const {
        return m_ints.m_int_2.phys_antisym_element(i, j, k, l);
    }

    defs::ham_t get_coeff(const ConnectionGeneral& conn) {
        switch (conn.get_signature()) {
            case ConnectionGeneral::signature(0, 0): return m_ints.m_int_0;
            case ConnectionGeneral::signature(1, 1):
                return m_ints.m_int_1.get(conn.m_fermion.m_cre[0], conn.m_fermion.m_ann[0]);
            case ConnectionGeneral::signature(2, 2):
                return m_ints.m_int_2.get(
                        conn.m_fermion.m_cre[0],
                        conn.m_fermion.m_cre[1],
                        conn.m_fermion.m_ann[0],
                        conn.m_fermion.m_ann[1]);
            default: return 0;
        }
    }
};

struct HubbardHamData : HamData {
    FermionIntegralArrays m_ints;

    HubbardHamData(FermionIntegralArrays &&ints) :
            HamData(ints.m_nsite,
                    conserve_fermion_exlvl(1, 1),
                    conserve_fermion_exlvl(0, 0),
                    null_exlvl(), null_exlvl(), null_exlvl()), m_ints(std::move(ints)) {}

    defs::ham_t get_coeff(const ConnectionGeneral& conn) {
        switch (conn.get_signature()) {
            case ConnectionGeneral::signature(0, 0): return m_ints.m_int_0;
            case ConnectionGeneral::signature(1, 1):
                return m_ints.m_int_1.get(conn.m_fermion.m_cre[0], conn.m_fermion.m_ann[0]);
            case ConnectionGeneral::signature(2, 2):
                return m_ints.m_int_2.get(
                        conn.m_fermion.m_cre[0],
                        conn.m_fermion.m_cre[1],
                        conn.m_fermion.m_ann[0],
                        conn.m_fermion.m_ann[1]);
            default: return 0;
        }
    }
};


struct HubbardHolstein : HamData {
    FermionIntegralArrays m_ints;
    const double m_v;
    const double m_omega;

    HubbardHolstein(FermionIntegralArrays &&ints, double v, double omega) :
            HamData(ints.m_nsite,
                    conserve_fermion_exlvl(1, 1),
                    conserve_fermion_exlvl(0, 0),
                    {{0, 0, 0, 1},
                     {0, 0, 0, 1}},
                    {{0, 0, 1, 0},
                     {0, 0, 1, 0}},

                    null_exlvl(), null_exlvl()), m_ints(std::move(ints)), m_v(v), m_omega(omega) {}

    defs::ham_t get_coeff_0000() {
        return m_ints.m_int_0;
    }

    defs::ham_t get_coeff_1100(const size_t &icre, const size_t &iann) {
        return m_ints.m_int_1.get(icre, iann);
    }

    defs::ham_t get_coeff_2200(const size_t &icre, const size_t &jcre, const size_t &iann, const size_t &jann) {
        return m_ints.m_int_2.phys_antisym_element(icre, jcre, iann, jann);
    }

    defs::ham_t get_coeff_1101(const size_t &icre_f, const size_t &iann_f, const size_t iann_b) {
        const size_t icre_site_f = icre_f < m_nsite ? icre_f : icre_f - m_nsite;
        const size_t iann_site_f = iann_f < m_nsite ? iann_f : iann_f - m_nsite;
        return (icre_site_f == iann_site_f && iann_b == icre_site_f) ? m_v : 0.0;
    }

    defs::ham_t get_coeff_1110(const size_t &icre_f, const size_t &iann_f, const size_t icre_b) {
        const size_t icre_site_f = icre_f < m_nsite ? icre_f : icre_f - m_nsite;
        const size_t iann_site_f = iann_f < m_nsite ? iann_f : iann_f - m_nsite;
        return (icre_site_f == iann_site_f && icre_b == icre_site_f) ? m_v : 0.0;
    }

    defs::ham_t get_coeff_0011(const size_t iann_b, const size_t icre_b) {
        return iann_b == icre_b ? m_omega : 0.0;
    }
};

/*
 * "coeff": single ham term from a promoted connection
 * "element": full matrix element between two ONVs
 */


template<typename data_t>
struct Ham {
    static_assert(std::is_base_of<HamData, data_t>::value, "Template arg must be derived from HamData");
    const data_t m_data;

    Ham(data_t &&data) : m_data(std::move(data)) {}

    const size_t &nsite() const {
        return static_cast<const HamData &>(m_data).m_nsite;
    }

    defs::ham_t get_element_0000() {

    }

    defs::ham_t get_coeff(const ConnectionGeneral& conn) {
        switch (conn.get_signature()) {
            case ConnectionGeneral::signature(0, 0):
                return m_data.get_coeff_0000();
            case ConnectionGeneral::signature(1, 1):
                return m_data.get_coeff_1100(conn.m_fermion.m_cre[0], conn.m_fermion.m_ann[0]);
            case ConnectionGeneral::signature(2, 2):
                return m_data.get_coeff_2200(
                        conn.m_fermion.m_cre[0],
                        conn.m_fermion.m_cre[1],
                        conn.m_fermion.m_ann[0],
                        conn.m_fermion.m_ann[1]);
            default: return 0;
        }
    }

    defs::ham_t get_element(const ConnectionGeneral& conn) {
        switch (conn.get_signature()) {
        }
    }

};
#endif


/*
 * hamiltonian has many HamTerms
 * given a Connection, the Hamiltonian matrix element can be computed
 *
 * get_element_0: (diagonal)
 * get_element_1: (single fermionic excitation)
 * get_element_1: (single excitation)
 */





/*
struct ConnectionGeneral {
    size_t m_cre_f;
    size_t m_ann_f;


};

struct HamiltonianGeneral {
    const std::vector<std::array<size_t, 4>> m_max_orders;
    // 4 operator types + 2 excit levels
    NdArray<bool, 6> m_nonzero_orders;
    const std::vector<defs::inds> m_exlvls;
    //HamiltonianGeneral(std::vector<std::array<size_t, 4>>& m_orders)
    HamiltonianGeneral(){
        m_nonzero_orders.m_format
    }

};
*/
#endif

#endif //M7_HAMILTONIAN_H
