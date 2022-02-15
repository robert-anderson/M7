//
// Created by anderson on 2/9/22.
//

#include "ForeachVirtual.h"

const foreach_virtual::rtnd::inds_t &foreach_virtual::rtnd::Base::value() const { return m_value; }

const size_t &foreach_virtual::rtnd::Base::iiter() const { return m_iiter; }

size_t foreach_virtual::rtnd::Base::ind_sum() const {
    return std::accumulate(m_value.cbegin(), m_value.cend(), 0ul);
}

foreach_virtual::rtnd::Base::Base(size_t nind, size_t nterm) : m_value(nind, 0), m_nind(nind), m_niter(nterm) {}

void foreach_virtual::rtnd::Base::loop() {
    if (m_nind) {
        m_iiter = ~0ul;
        try {
            throwing_loop();
            DEBUG_ASSERT_EQ(m_iiter+1, m_niter, "loop completed after incorrect number of iterations");
        }
        catch (const ExitLoop &) {}
    } else m_iiter = 0ul;
}

size_t foreach_virtual::rtnd::Unrestricted::nterm(const foreach_virtual::rtnd::inds_t &shape) {
    if (shape.empty()) return 0;
    size_t n = 1;
    for (const auto &i: shape) n *= i;
    return n;
}

void foreach_virtual::rtnd::Unrestricted::level_loop(size_t ilevel) {
    const auto &iind = ilevel - 1;
    auto &ind = m_value[iind];
    const auto &extent = m_shape[iind];
    try {
        if (ilevel < m_nind)
            for (ind = 0ul; ind < extent; ++ind) level_loop(ilevel + 1);
        else {
            for (ind = 0ul; ind < extent; ++ind) {
                ++m_iiter;
                body();
            }
        }
    }
    catch (const ExitLoop& ex){throw ex;}
}

foreach_virtual::rtnd::Unrestricted::Unrestricted(const foreach_virtual::rtnd::inds_t &shape) :
    Base(shape.size(), nterm(shape)), m_shape(shape) {}

foreach_virtual::rtnd::Unrestricted::Unrestricted(const size_t &nind, const size_t &extent) :
    Unrestricted(std::vector<size_t>(nind, extent)) {}

void foreach_virtual::rtnd::Unrestricted::throwing_loop() {
    m_iiter = ~0ul;
    level_loop(1);
}
