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
    const size_t m_size;
    /**
     * number of physical dimensions
     */
    const size_t m_ndim;
    const std::vector<std::vector<int>> m_momvecs;
    mutable std::vector<int> m_momvec_work;
    mutable defs::inds m_inds_work;


    static defs::inds make_momentum_shape(const defs::inds& wave_shape);
    static std::vector<std::vector<int>> make_momvecs(const defs::inds& wave_shape);

public:
    static size_t size(const defs::inds& wave_shape);
    static size_t size(size_t ndim, size_t nwave);

    explicit Planewaves(const defs::inds& wave_shape);

    Planewaves(size_t ndim, size_t nwave);

    const std::vector<int>& operator[](const size_t& i) const;

    bool conserving(const size_t& icre, const size_t& iann) const;

    bool conserving(const size_t& icre1, const size_t& icre2, const size_t& iann1, const size_t& iann2) const;

    int diff(const size_t& idim, const size_t& i, const size_t& j) const;

    void diff(std::vector<int>& diff, const size_t& i, const size_t& j) const;

    size_t diff(const size_t& i, const size_t& j) const;

    size_t encode(const std::vector<int>& momvec) const;

    const size_t& size() const;

    int kinetic_energy(const size_t& i) const;
};


#endif //M7_PLANEWAVES_H
