//
// Created by rja on 12/12/22.
//

#ifndef M7_INSTESTS_H
#define M7_INSTESTS_H

#include "CommutingObservable.h"
#include "M7_lib/wavefunction/Reference.h"

/**
 * Instantaneous estimators. Optional stats that are linear in the wavefunction and can therefore be timestep-resolved
 */
struct InstEsts {
    std::unique_ptr<commuting_obs::SpinSquare> m_spin_square = nullptr;
    InstEsts(const sys::Sector sector, const wf::References* refs, const conf::InstEsts& opts);

    void begin_cycle(uint_t icycle);

    void make_numerator_contribs(const Walker& walker);

    void end_cycle(uint_t icycle);
};


#endif //M7_INSTESTS_H
