//
// Created by Robert J. Anderson on 08/03/2022.
//

#include "Planewaves.h"

defs::uintv_t Planewaves::make_momentum_shape(const defs::uintv_t &wave_shape) {
    defs::uintv_t shape;
    for (auto extent: wave_shape) shape.push_back(2*extent+1);
    return shape;
}
std::vector<std::vector<int>> Planewaves::make_momvecs(const defs::uintv_t &wave_shape) {
    std::vector<std::vector<int>> momvecs;
    auto momentum_shape = make_momentum_shape(wave_shape);
    momvecs.reserve(size(wave_shape));
    auto fn = [&wave_shape, &momvecs](const defs::uintv_t& inds){
        momvecs.emplace_back();
        for (size_t idim=0ul; idim < wave_shape.size(); ++idim){
            momvecs.back().push_back(int(inds[idim]) - wave_shape[idim]);
        }
    };
    basic_foreach::rtnd::Unrestricted foreach(momentum_shape);
    foreach.loop(fn);
    DEBUG_ASSERT_EQ(momvecs.size(), momvecs.capacity(), "incorrect number of momentum vectors generated");
    return momvecs;
}

size_t Planewaves::size(const defs::uintv_t &wave_shape) {
    const auto shape = make_momentum_shape(wave_shape);
    return NdFormatD(shape).m_nelement;
}

size_t Planewaves::size(size_t ndim, size_t nwave) {
    return size(defs::uintv_t(ndim, nwave));
}

Planewaves::Planewaves(const defs::uintv_t& wave_shape) :
    m_wave_format(wave_shape), m_momentum_format(make_momentum_shape(wave_shape)),
    m_size(m_momentum_format.m_nelement), m_ndim(m_momentum_format.m_nind),
    m_momvecs(make_momvecs(wave_shape)){}

Planewaves::Planewaves(size_t ndim, size_t nwave) : Planewaves(defs::uintv_t(ndim, nwave)){}

const std::vector<int> &Planewaves::operator[](const size_t &i) const {
    DEBUG_ASSERT_LT(i, m_size, "basis function index OOB");
    return m_momvecs[i];
}

bool Planewaves::conserving(const size_t &icre, const size_t &iann) const {
    auto& cre_mom = (*this)[icre];
    auto& ann_mom = (*this)[iann];
    for (size_t idim=0ul; idim<m_ndim; ++idim) if (cre_mom[idim]!=ann_mom[idim]) return false;
    return true;
}

bool Planewaves::conserving(const size_t &icre1, const size_t &icre2, const size_t &iann1, const size_t &iann2) const {
    auto& cre_mom1 = (*this)[icre1];
    auto& cre_mom2 = (*this)[icre2];
    auto& ann_mom1 = (*this)[iann1];
    auto& ann_mom2 = (*this)[iann2];
    for (size_t idim=0ul; idim<m_ndim; ++idim){
        if (cre_mom1[idim]+cre_mom2[idim]!=ann_mom1[idim]+ann_mom2[idim]) return false;
    }
    return true;
}

int Planewaves::diff(const size_t &idim, const size_t &i, const size_t &j) const {
    DEBUG_ASSERT_LT(idim, m_ndim, "dimension index OOB");
    auto& imom = (*this)[i];
    auto& jmom = (*this)[j];
    return imom[idim]-jmom[idim];
}

void Planewaves::diff(std::vector<int> &diff, const size_t &i, const size_t &j) const {
    auto& imom = (*this)[i];
    auto& jmom = (*this)[j];
    diff.clear();
    for (size_t idim=0ul; idim<m_ndim; ++idim) diff.push_back(imom[idim]-jmom[idim]);
}

size_t Planewaves::diff(const size_t &i, const size_t &j) const {
    diff(m_momvec_work, i, j);
    return encode(m_momvec_work);
}

size_t Planewaves::encode(const std::vector<int> &momvec) const {
    m_inds_work.clear();
    for (size_t idim=0ul; idim<m_ndim; ++idim)
        m_inds_work.push_back(momvec[idim]+m_wave_format.m_shape[idim]);
    return m_momentum_format.flatten(m_inds_work);
}

const size_t &Planewaves::size() const {
    return m_size;
}

int Planewaves::kinetic_energy(const size_t &i) const {
    auto& momvec = (*this)[i];
    int tot = 0;
    for (auto& mom: momvec) tot+=mom*mom;
    return tot;
}
