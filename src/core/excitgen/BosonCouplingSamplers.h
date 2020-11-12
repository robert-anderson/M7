//
// Created by rja on 12/11/2020.
//

#ifndef M7_BOSONCOUPLINGSAMPLERS_H
#define M7_BOSONCOUPLINGSAMPLERS_H

#include <src/core/basis/DecodedDeterminant.h>
#include <src/core/basis/Connections.h>
#include "src/core/field/Views.h"

struct BosonCouplingSamplers {

    size_t m_nboson_max;
    BosonCouplingSamplers(size_t nboson_max):m_nboson_max(nboson_max){}

    bool draw_single(const views::FermiBosOnv &src_onv, views::FermiBosOnv &dst_onv,
                     const OccupiedOrbitals &occ, defs::prob_t &prob, defs::ham_t &helem, conn::AsFermiBosOnv &anticonn) {
        /*
         * TODO: James
         *  make sure to either prevent the draw of a boson onv which exceeds the cutoff, or at
         *  least return false in such cases
         */
        return false;
    }

};


#endif //M7_BOSONCOUPLINGSAMPLERS_H
