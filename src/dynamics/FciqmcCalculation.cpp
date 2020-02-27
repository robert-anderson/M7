//
// Created by Robert John Anderson on 2020-02-11.
//

#include "FciqmcCalculation.h"
#include "src/io/Logging.h"

void FciqmcCalculation::execute(size_t ncycle){
    log::write("Starting FCIQMC main loop.");
    for(size_t icycle=0ul; icycle<ncycle; ++icycle){
        m_psi.propagate(m_prop);
        m_psi.communicate();
        m_psi.consolidate_incoming_weight();
        m_psi.annihilate();
        m_prop.update(icycle, m_psi.norm_growth_rate());
        m_psi.write_iter_stats();
    }
}
