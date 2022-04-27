//
// Created by Robert J. Anderson on 04/03/2022.
//

#ifndef M7_DECODEDFRMBOSONV_H
#define M7_DECODEDFRMBOSONV_H

#include "DecodedMbf.h"

namespace decoded_mbf {
    namespace frmbos {

        struct Base : Cache {
        protected:
            const FrmBosOnvField &m_mbf;
        public:
            Base(const FrmBosOnvField &mbf);

            bool is_valid() const override;
        };

        /**
         * used for either spinorb, site, or mode indices
         */
        struct SimpleBase : Base, SimpleContainer {
        protected:
            const defs::inds &validated() const;

        public:
            explicit SimpleBase(const FrmBosOnvField &mbf);
        };

        /**
         * site indices with any electrons and any associated bosons (imode==isite) e.g. holstein
         */
        struct OccSitesNonzeroBosons : SimpleBase {
            explicit OccSitesNonzeroBosons(const FrmBosOnvField &mbf);

            const defs::inds &get();

        };
    }


    struct FrmBosOnv {
        const FrmBosOnvField &m_mbf;
        frmbos::OccSitesNonzeroBosons m_occ_sites_nonzero_bosons;

        FrmBosOnv(const FrmBosOnvField &mbf);

        /**
         * clear all cached assets
         */
        void clear();
    };
}

#endif //M7_DECODEDFRMBOSONV_H
