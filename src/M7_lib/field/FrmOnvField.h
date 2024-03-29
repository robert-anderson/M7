//
// Created by Robert J. Anderson on 07/04/2021.
//

#ifndef M7_FRMONVFIELD_H
#define M7_FRMONVFIELD_H

#include <M7_lib/basis/BasisData.h>
#include <M7_lib/caches/DecodedFrmOnv.h>

#include "BitsetField.h"

struct FrmOnvField : BitsetField<uint_t, 2> {
    typedef BitsetField<uint_t, 2> base_t;
    using base_t::get;
    using base_t::set;
    using base_t::clr;
    using base_t::put;
    using base_t::operator=;
    using base_t::inds_t;

    const sys::frm::Basis m_basis;
    /**
     * a refreshable cache of useful representations for excitation generation and enumeration
     */
    mutable decoded_mbf::FrmOnv m_decoded;

private:
    /**
     * number of data words required for the storage of a single spin channel (alpha, beta)
     */
    const uint_t m_dsize_spin_channel;

public:

    FrmOnvField(Row* row, const sys::frm::Basis& basis, str_t name="");
    /*
     * all ONVs implement the following ctor
     */
    FrmOnvField(Row* row, const sys::Basis& basis, str_t name="");
    /*
     * this particular MBF only needs the basis, but future MBF types might need the full sector information, and so
     * a common interface is realised by implementing a ctor of the following form in all MBFs
     */
    FrmOnvField(Row* row, const sys::frm::Sector& sector, str_t name="");
    FrmOnvField(Row* row, const sys::Sector& sector, str_t name="");
    FrmOnvField(Row* row, const FrmOnvField& other);
    FrmOnvField(const FrmOnvField& other);

    FrmOnvField& operator=(const FrmOnvField& other) {
        base_t::operator=(other);
        return *this;
    }

    FrmOnvField &operator=(std::pair<const uintv_t &, const uintv_t &> setbits);

    bool operator==(const FrmOnvField& other) const {
        return base_t::operator==(other);
    }

    bool operator==(const uintv_t& inds) const;

    void set(const uint_t& bit_offset, const uintv_t& setbits);

    void set(const uint_t& site_offset, const uintv_t& setbits_alpha, const uintv_t& setbits_beta);

    void set(const uintv_t& setbits_alpha, const uintv_t& setbits_beta);

    void set_spins(const uintv_t& alpha_sites);

    void put_spin_channel(const uint_t& ispin, bool set);

    void clr(const uint_t& bit_offset, const uintv_t& clrbits);

    void clr(const uint_t& site_offset, const uintv_t& clrbits_alpha, const uintv_t& clrbits_beta);

    void excite(const uint_t &i, const uint_t &j) {
        clr(tbegin(), i);
        set(tbegin(), j);
    }
    void excite(inds_t ann, inds_t cre){
        clr(tbegin(), ann);
        set(tbegin(), cre);
    }

    void excite(const uint_t &i, const uint_t &j, const uint_t &k, const uint_t &l) {
        clr(tbegin(), i);
        clr(tbegin(), j);
        set(tbegin(), k);
        set(tbegin(), l);
    }

    void excite(inds_t ann1, inds_t ann2, inds_t cre1, inds_t cre2) {
        clr(tbegin(), ann1);
        clr(tbegin(), ann2);
        set(tbegin(), cre1);
        set(tbegin(), cre2);
    }

    /**
     * 2 x the total z-axis projection of total spin
     * @return
     *  nalpha-nbeta where alpha is conventionally the 0 spin channel
     */
    int ms2() const;

    uint_t site_nocc(const uint_t& isite) const;

    bool all_sites_single_occ() const;

    str_t to_string() const override;

    /**
     * @param idataword
     *  dataword index within the spin channel
     * @return
     *  raw dataword representing a portion of the alpha spin channel of the FrmOnv
     */
    uint_t get_alpha_dataword(uint_t idataword) const;

    /**
     * slightly more involved than the alpha channel counterpart. since the spin channels are not separated by a word
     * boundary, the beta channel dataword must in general be formed by shift-AND operations on adjacent datawords of
     * the buffer.
     * @param idataword
     *  dataword index within the spin channel
     * @return
     *  raw dataword representing a portion of the beta spin channel of the FrmOnv
     */
    uint_t get_beta_dataword(uint_t idataword) const;

    template<typename body_fn_t>
    void foreach_alpha(const body_fn_t& fn) const {
        auto get_work_fn = [&](uint_t idataword){
            return get_alpha_dataword(idataword);
        };
        setbit_foreach::single<uint_t>(m_dsize_spin_channel, fn, get_work_fn);
    }
    uint_t nalpha() const;

    template<typename body_fn_t>
    void foreach_beta(const body_fn_t& fn) const {
        auto get_work_fn = [&](uint_t idataword){
            return get_beta_dataword(idataword);
        };
        setbit_foreach::single<uint_t>(m_dsize_spin_channel, fn, get_work_fn);
    }
    uint_t nbeta() const;


    /**
     * efficiently iterate over the singly-occupied site indices
     * @tparam body_fn_t
     *  functor to call each time a singly-occupied site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_open_shell(const body_fn_t& fn) const {
        auto get_work_fn = [&](uint_t idataword){
            return get_alpha_dataword(idataword) ^ get_beta_dataword(idataword);
        };
        setbit_foreach::single<uint_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    uint_t nopen_shell() const;

    /**
     * efficiently iterate over the indices of those sites where the alpha spin orbital is occupied but the beta spin
     * orbital is vacant
     * @tparam body_fn_t
     *  functor to call each time an alpha-occupied, beta-vacant site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_open_shell_alpha(const body_fn_t& fn) const {
        auto get_work_fn = [&](uint_t idataword){
            return get_alpha_dataword(idataword) &~ get_beta_dataword(idataword);
        };
        setbit_foreach::single<uint_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    uint_t nopen_shell_alpha() const;

    /**
     * efficiently iterate over the indices of those sites where the alpha spin vacant is occupied but the beta spin
     * orbital is occupied
     * @tparam body_fn_t
     *  functor to call each time an alpha-occupied, beta-vacant site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_open_shell_beta(const body_fn_t& fn) const {
        auto get_work_fn = [&](uint_t idataword){
            return ~get_alpha_dataword(idataword) & get_beta_dataword(idataword);
        };
        setbit_foreach::single<uint_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    uint_t nopen_shell_beta() const;

    /**
     * efficiently iterate over the doubly-occupied site indices
     * @tparam body_fn_t
     *  functor to call each time a doubly-occupied site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_closed_shell(const body_fn_t& fn) const {
        auto get_work_fn = [&](uint_t idataword){
            return get_alpha_dataword(idataword) & get_beta_dataword(idataword);
        };
        setbit_foreach::single<uint_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    uint_t nclosed_shell() const;

    /**
     * efficiently iterate over the site indices with any non-zero occupation
     * @tparam body_fn_t
     *  functor to call each time an occupied site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_occupied_site(const body_fn_t& fn) const {
        auto get_work_fn = [&](uint_t idataword){
            return get_alpha_dataword(idataword) | get_beta_dataword(idataword);
        };
        setbit_foreach::single<uint_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    uint_t noccupied_site() const;

    /**
     * efficiently iterate over the unoccupied site indices
     * @tparam body_fn_t
     *  functor to call each time a unoccupied site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_unoccupied_site(const body_fn_t& fn) const {
        uint_t isite = 0ul;
        auto tmp = [&isite, &fn](uint_t iocc_site){
            for (; isite<iocc_site; ++isite){
                fn(isite);
            }
            // step over the occupied site
            ++isite;
        };
        foreach_occupied_site(tmp);
        // the remaining sites are also unoccupied
        for (;isite<m_basis.m_nsite; ++isite){
            fn(isite);
        }
    }

    uint_t nunoccupied_site() const;

};


#endif //M7_FRMONVFIELD_H
