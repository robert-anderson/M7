//
// Created by Robert J. Anderson on 12/8/21.
//

#ifndef M7_GENERALFRMHAM_H
#define M7_GENERALFRMHAM_H

#include "FrmHam.h"
#include "M7_lib/integrals/IntegralArray1e.h"
#include "M7_lib/integrals/IntegralArray2e.h"
#include "M7_lib/io/FcidumpTextFileReader.h"


struct IntegralReader {
    struct IterData {
        uintv_t m_inds = uintv_t(4ul, ~0ul);
        ham_t m_value {};
        uint_t m_ranksig {};
        uint_t m_exsig {};
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

struct Hdf5IntegralReader : IntegralReader {
    hdf5::FileReader m_reader;
    Hdf5IntegralReader(const FcidumpInfo& info, bool /*spin_major*/): m_reader(info.m_fname){
        std::vector<ham_t> values_2e;
        std::cout << m_reader.child_name(0) << std::endl;
        std::cout << m_reader.child_name(1) << std::endl;
        std::cout << m_reader.child_name(2) << std::endl;
        std::cout << m_reader.child_name(3) << std::endl;
        std::cout << m_reader.child_name(4) << std::endl;
        std::cout << m_reader.child_name(5) << std::endl;

        REQUIRE_TRUE(m_reader.child_exists("TWO_EL_INT_VALUES"), "2e integral values not found in HDF5 file");
        m_reader.load("TWO_EL_INT_VALUES", values_2e);
        std::cout << values_2e << std::endl;
    }
    bool next(IterData& /*data*/) override {
        return false;
    }

    void goto_first_1e() override {
    }

    void goto_first_2e() override {
    }

    ham_t ecore() const override {
        return m_reader.read_attr<double>("CORE_ENERGY");
    }

    bool complex_valued() const override {
        return false;
    }


};

//struct Hdf5IntegralReader : IntegralReader {
//    hdf5::FileReader m_reader;
//
//    Hdf5IntegralReader(const std::string& fname, bool spin_major): m_reader(fname){}
//
//    bool next(IterData& data) override;
//
//    void goto_first_1e() override;
//
//    void goto_first_2e() override;
//
//    ham_t ecore() const override;
//
//    bool complex_valued() const override;
//};

struct GeneralFrmHam : FrmHam {

    typedef integrals_1e::Array<ham_t> ints_1e_t;
    typedef integrals_2e::Array<ham_t> ints_2e_t;
    typedef std::unique_ptr<ints_1e_t> ints_1e_ptr_t;
    typedef std::unique_ptr<ints_2e_t> ints_2e_ptr_t;

    struct Integrals {
        ints_1e_ptr_t m_1e;
        ints_2e_ptr_t m_2e;
    };

private:
    const FcidumpInfo m_info;

    static void log_ints_sym(integrals_1e::syms::Sym sym, bool initial);

    static void log_ints_sym(integrals_2e::syms::Sym sym, bool initial);

    Integrals make_ints(IntegralReader& reader);

    Integrals make_ints(const FcidumpInfo& info, bool spin_major) {
        const std::string fmt = "Reading fermion Hamiltonian coefficients from {} file \"" + info.m_fname + "\"...";
        if (info.m_impl==FcidumpInfo::CSV) {
            log::info(fmt, "plain text CSV");
            CsvIntegralReader reader(info, spin_major);
            return make_ints(reader);
        }
        else if (info.m_impl==FcidumpInfo::HDF5) {
            log::info(fmt, "HDF5 binary");
            Hdf5IntegralReader reader(info, spin_major);
            return make_ints(reader);
        }
        return {nullptr, nullptr};
    }

public:

    const Integrals m_ints;

    GeneralFrmHam(const FcidumpInfo& info, bool spin_major);

    explicit GeneralFrmHam(opt_pair_t opts);

    ham_t get_coeff_1100(uint_t a, uint_t i) const override;

    ham_t get_coeff_2200(uint_t a, uint_t b, uint_t i, uint_t j) const override;


    ham_t get_element_0000(const field::FrmOnv &onv) const override;

    ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    /**
     * @param prng
     *  random number generator
     * @param opts
     *  propagator options (so excitation generation settings can be read from the configuration document)
     * @param h
     *  Fermion H object to be **treated** as a GeneralFrmHam
     * @return
     *  list of excit gens required to stochastically propagate the given h
     */
    static excit_gen_list_t make_excit_gens(PRNG& prng, const conf::Propagator& opts, const FrmHam& h);

    excit_gen_list_t make_excit_gens(PRNG &prng, const conf::Propagator& opts) const override;

    conn_foreach_list_t make_foreach_iters() const override;

    uint_t default_nelec() const override;

    int default_ms2_value() const override;
};


#endif //M7_GENERALFRMHAM_H
