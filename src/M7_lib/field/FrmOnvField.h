//
// Created by Robert J. Anderson on 07/04/2021.
//

#ifndef M7_FRMONVFIELD_H
#define M7_FRMONVFIELD_H

#include <M7_lib/basis/BasisData.h>
#include <M7_lib/caches/DecodedFrmOnv.h>

#include "BitsetField.h"

struct FrmOnvField : BitsetField<size_t, 2> {
    typedef BitsetField<size_t, 2> base_t;
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
    const size_t m_dsize_spin_channel;

public:

    FrmOnvField(Row* row, const sys::frm::Basis& basis, std::string name="");
    /*
     * all ONVs implement the following ctor
     */
    FrmOnvField(Row* row, const sys::Basis& basis, std::string name="");
    /*
     * this particular MBF only needs the basis, but future MBF types might need the full sector information, and so
     * a common interface is realised by implementing a ctor of the following form in all MBFs
     */
    FrmOnvField(Row* row, const sys::frm::Sector& sector, std::string name="");
    FrmOnvField(Row* row, const sys::Sector& sector, std::string name="");
    FrmOnvField(const FrmOnvField& other);

    FrmOnvField& operator=(const FrmOnvField& other) {
        base_t::operator=(other);
        return *this;
    }

    FrmOnvField &operator=(std::pair<const defs::inds_t &, const defs::inds_t &> setbits);

    bool operator==(const FrmOnvField& other) const {
        return base_t::operator==(other);
    }

    bool operator==(const defs::inds_t& inds) const;

    void set(const size_t& bit_offset, const defs::inds_t& setbits);

    void set(const size_t& site_offset, const defs::inds_t& setbits_alpha, const defs::inds_t& setbits_beta);

    void set(const defs::inds_t& setbits_alpha, const defs::inds_t& setbits_beta);

    void set_spins(const defs::inds_t& alpha_sites);

    void put_spin_channel(const size_t& ispin, bool set);

    void clr(const size_t& bit_offset, const defs::inds_t& clrbits);

    void clr(const size_t& site_offset, const defs::inds_t& clrbits_alpha, const defs::inds_t& clrbits_beta);

    void excite(const size_t &i, const size_t &j) {
        auto* dptr = reinterpret_cast<size_t *>(begin());
        clr(dptr, i);
        set(dptr, j);
    }
    void excite(inds_t ann, inds_t cre){
        auto* dptr = reinterpret_cast<size_t *>(begin());
        clr(dptr, ann);
        set(dptr, cre);
    }

    void excite(const size_t &i, const size_t &j, const size_t &k, const size_t &l) {
        auto* dptr = reinterpret_cast<size_t *>(begin());
        clr(dptr, i);
        clr(dptr, j);
        set(dptr, k);
        set(dptr, l);
    }
    void excite(inds_t ann1, inds_t ann2, inds_t cre1, inds_t cre2) {
        auto* dptr = reinterpret_cast<size_t *>(begin());
        clr(dptr, ann1);
        clr(dptr, ann1);
        set(dptr, cre1);
        set(dptr, cre2);
    }

    /**
     * 2 x the total z-axis projection of total spin
     * @return
     *  nalpha-nbeta where alpha is conventionally the 0 spin channel
     */
    int ms2() const;

    size_t site_nocc(const size_t& isite) const;

    bool all_sites_single_occ() const;

    std::string to_string() const override;

    /**
     * @param idataword
     *  dataword index within the spin channel
     * @return
     *  raw dataword representing a portion of the alpha spin channel of the FrmOnv
     */
    size_t get_alpha_dataword(size_t idataword) const;

    /**
     * slightly more involved than the alpha channel counterpart. since the spin channels are not separated by a word
     * boundary, the beta channel dataword must in general be formed by shift-AND operations on adjacent datawords of
     * the buffer.
     * @param idataword
     *  dataword index within the spin channel
     * @return
     *  raw dataword representing a portion of the beta spin channel of the FrmOnv
     */
    size_t get_beta_dataword(size_t idataword) const;

    template<typename body_fn_t>
    void foreach_alpha(const body_fn_t& fn) const {
        auto get_work_fn = [&](size_t idataword){
            return get_alpha_dataword(idataword);
        };
        setbit_foreach::single<size_t>(m_dsize_spin_channel, fn, get_work_fn);
    }
    size_t nalpha() const;

    template<typename body_fn_t>
    void foreach_beta(const body_fn_t& fn) const {
        auto get_work_fn = [&](size_t idataword){
            return get_beta_dataword(idataword);
        };
        setbit_foreach::single<size_t>(m_dsize_spin_channel, fn, get_work_fn);
    }
    size_t nbeta() const;


    /**
     * efficiently iterate over the singly-occupied site indices
     * @tparam body_fn_t
     *  functor to call each time a singly-occupied site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_open_shell(const body_fn_t& fn) const {
        auto get_work_fn = [&](size_t idataword){
            return get_alpha_dataword(idataword) ^ get_beta_dataword(idataword);
        };
        setbit_foreach::single<size_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    size_t nopen_shell() const;

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
        auto get_work_fn = [&](size_t idataword){
            return get_alpha_dataword(idataword) &~ get_beta_dataword(idataword);
        };
        setbit_foreach::single<size_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    size_t nopen_shell_alpha() const;

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
        auto get_work_fn = [&](size_t idataword){
            return ~get_alpha_dataword(idataword) & get_beta_dataword(idataword);
        };
        setbit_foreach::single<size_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    size_t nopen_shell_beta() const;

    /**
     * efficiently iterate over the doubly-occupied site indices
     * @tparam body_fn_t
     *  functor to call each time a doubly-occupied site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_closed_shell(const body_fn_t& fn) const {
        auto get_work_fn = [&](size_t idataword){
            return get_alpha_dataword(idataword) & get_beta_dataword(idataword);
        };
        setbit_foreach::single<size_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    size_t nclosed_shell() const;

    /**
     * efficiently iterate over the site indices with any non-zero occupation
     * @tparam body_fn_t
     *  functor to call each time an occupied site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_occupied_site(const body_fn_t& fn) const {
        auto get_work_fn = [&](size_t idataword){
            return get_alpha_dataword(idataword) | get_beta_dataword(idataword);
        };
        setbit_foreach::single<size_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    size_t noccupied_site() const;

    /**
     * efficiently iterate over the unoccupied site indices
     * @tparam body_fn_t
     *  functor to call each time a unoccupied site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_unoccupied_site(const body_fn_t& fn) const {
        size_t isite = 0ul;
        auto tmp = [&isite, &fn](size_t iocc_site){
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

    size_t nunoccupied_site() const;

};


#endif //M7_FRMONVFIELD_H
