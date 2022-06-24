//
// Created by Robert J. Anderson on 04/03/2022.
//

#ifndef M7_DECODEDBOSONV_H
#define M7_DECODEDBOSONV_H

#include "DecodedMbf.h"

namespace decoded_mbf {
    namespace bos {

        struct Base : Cache {
        protected:
            const BosOnvField &m_mbf;
        public:
            Base(const BosOnvField &mbf);

            bool is_valid() const override;
        };

        struct SimpleBase : Base, SimpleContainer {
        protected:
            const defs::uintv_t &validated() const;

        public:
            explicit SimpleBase(const BosOnvField &mbf);
        };

        /**
         * boson mode indices with repetition
         *  e.g. [0, 2, 0, 3, 1] decodes as:
         *      [1, 1, 3, 3, 3, 4]
         */
        struct Expanded : SimpleBase {
            Expanded(const BosOnvField &mbf);

            const defs::uintv_t &get();
        };

        /**
         * indices of boson modes with any non-zero occupation
        *  e.g. [0, 2, 0, 3, 1] decodes as:
        *      [1, 3, 4]
         */
        struct OccModes : SimpleBase {
            OccModes(const BosOnvField &mbf);

            const defs::uintv_t &get();
        };
    }

    struct BosOnv {
        /**
         * all boson operators in product with repetitions
         */
        bos::Expanded m_expanded;
        /**
         * indices of all modes with non-zero occupation
         */
        bos::OccModes m_occ_modes;

        BosOnv(const BosOnvField &mbf);

        /**
         * clear all cached assets
         */
        void clear();
    };
}

#endif //M7_DECODEDBOSONV_H