//
// Created by Robert J. Anderson on 08/03/2022.
//

#ifndef M7_PLANEWAVES_H
#define M7_PLANEWAVES_H

#include <M7_lib/nd/NdFormatD.h>

/**
 * generic to fermionic and bosonic degrees of freedom
 * waves refers to the number of waves in each dimension
 * momenta refers to the possible values of wave vector
 * e.g. if number of waves is 3, then the number of momenta is 7: (-3, -2, -1, 0, 1, 2, 3)
 */
class Planewaves {
    const NdFormatD m_wave_format;
    const NdFormatD m_momentum_format;
    /**
     * number of degrees of freedom (sites or modes)
     */
    const uint_t m_size;
    /**
     * number of physical dimensions
     */
    const uint_t m_ndim;
    const v_t<v_t<int>> m_momvecs;
    mutable v_t<int> m_momvec_work;
    mutable uintv_t m_inds_work;


    static uintv_t make_momentum_shape(const uintv_t& wave_shape);
    static v_t<v_t<int>> make_momvecs(const uintv_t& wave_shape);

public:
    static uint_t size(const uintv_t& wave_shape);
    static uint_t size(uint_t ndim, uint_t nwave);

    explicit Planewaves(const uintv_t& wave_shape);

    Planewaves(uint_t ndim, uint_t nwave);

    const v_t<int>& operator[](const uint_t& i) const;

    bool conserving(const uint_t& icre, const uint_t& iann) const;

    bool conserving(const uint_t& icre1, const uint_t& icre2, const uint_t& iann1, const uint_t& iann2) const;

    int diff(const uint_t& idim, const uint_t& i, const uint_t& j) const;

    void diff(v_t<int>& diff, const uint_t& i, const uint_t& j) const;

    uint_t diff(const uint_t& i, const uint_t& j) const;

    uint_t encode(const v_t<int>& momvec) const;

    const uint_t& size() const;

    int kinetic_energy(const uint_t& i) const;
};


#endif //M7_PLANEWAVES_H
