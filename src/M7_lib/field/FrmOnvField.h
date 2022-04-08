//
// Created by rja on 07/04/2021.
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

    const FrmBasisData m_bd;
    /**
     * a refreshable cache of useful representations for excitation generation and enumeration
     */
    mutable decoded_mbf::FrmOnv m_decoded;

private:
    /**
     * number of data words required for the storage of a single spin channel (alpha, beta)
     */
    const size_t m_dsize_spin_channel;
    /**
     * alpha and beta strings in general share a dataword, so store the bit position which marks the interface
     */
    const size_t m_nbit_in_last_alpha_dataword;

public:

    FrmOnvField(Row* row, const FrmBasisData& bd, std::string name="");
    FrmOnvField(Row* row, size_t nsite, std::string name="");
    FrmOnvField(Row* row, const BasisData& bd, std::string name="");
    FrmOnvField(const FrmOnvField& other);

    FrmOnvField& operator=(const FrmOnvField& other) {
        base_t::operator=(other);
        return *this;
    }

    FrmOnvField &operator=(std::pair<const defs::inds &, const defs::inds &> setbits);

    bool operator==(const FrmOnvField& other) const {
        return base_t::operator==(other);
    }

    bool operator==(const defs::inds& inds) const;

    void set(const size_t& bit_offset, const defs::inds& setbits);

    void set(const size_t& site_offset, const defs::inds& setbits_alpha, const defs::inds& setbits_beta);

    void set(const defs::inds& setbits_alpha, const defs::inds& setbits_beta);

    void set_spins(const defs::inds& alpha_sites);

    void put_spin_channel(const size_t& ispin, bool set);

    void clr(const size_t& bit_offset, const defs::inds& clrbits);

    void clr(const size_t& site_offset, const defs::inds& clrbits_alpha, const defs::inds& clrbits_beta);

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

    int nalpha() const;

    size_t site_nocc(const size_t& isite) const;

    bool all_sites_single_occ() const;

    std::string to_string() const override;

    static size_t isite(const size_t& ibit, const size_t& nsite);

    static size_t ispin(const size_t& ibit, const size_t& nsite);

    static size_t ibit(const size_t& ispin, const size_t& isite, const size_t& nsite);

    static size_t ibit(std::pair<size_t, size_t> pair, const size_t& nsite);;

    /**
     * convenient conversion between spin-site and site indices
     * @param ibit
     *  bit index
     * @return
     *  associated site index
     */
    size_t isite(const size_t& ibit) const;
    /**
     * convenient conversion between spin-site and spin indices
     * @param ibit
     *  bit index
     * @return
     *  associated spin index (0 or 1)
     */
    size_t ispin(const size_t& ibit) const;
    /**
     * convenient conversion between spin and site indices and the flat "bit" index
     * @param ispin
     *  spin channel index
     * @param isite
     *  site index
     * @return
     *  flat index od spin-site
     */
    size_t ibit(const size_t& ispin, const size_t& isite) const;

    size_t ibit(std::pair<size_t, size_t>& pair) const;

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

    template<typename body_fn_t>
    void foreach_beta(const body_fn_t& fn) const {
        auto get_work_fn = [&](size_t idataword){
            return get_beta_dataword(idataword);
        };
        setbit_foreach::single<size_t>(m_dsize_spin_channel, fn, get_work_fn);
    }


    /**
     * efficiently iterate over the singly-occupied site indices
     * @tparam body_fn_t
     *  functor to call each time a singly-occupied site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_openshell(const body_fn_t& fn) const {
        auto get_work_fn = [&](size_t idataword){
            return get_alpha_dataword(idataword) ^ get_beta_dataword(idataword);
        };
        setbit_foreach::single<size_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    /**
     * efficiently iterate over the indices of those sites where the alpha spin orbital is occupied but the beta spin
     * orbital is vacant
     * @tparam body_fn_t
     *  functor to call each time an alpha-occupied, beta-vacant site is found
     * @param fn
     *  body_fn_t instance
     */
    template<typename body_fn_t>
    void foreach_openshell_alpha(const body_fn_t& fn) const {
        auto get_work_fn = [&](size_t idataword){
            return get_alpha_dataword(idataword) &~ get_beta_dataword(idataword);
        };
        setbit_foreach::single<size_t>(m_dsize_spin_channel, fn, get_work_fn);
    }

    size_t nopenshell() const;

    size_t nopenshell_alpha() const;
};


#endif //M7_FRMONVFIELD_H
