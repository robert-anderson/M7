//
// Created by rja on 24/08/2021.
//

#ifndef M7_DECODEDMBF_H
#define M7_DECODEDMBF_H

#include <algorithm>
#include <M7_lib/basis/AbelianGroup.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/nd/NdFormat.h>
#include <M7_lib/defs.h>

struct FrmOnvField;
struct BosOnvField;
struct FrmBosOnvField;

/**
 * provides a store which only updates the flat and spin/sym-partitioned occupied and vacant spin orbital indices when
 * needed
 */
namespace decoded_mbf {

    /**
     * base class for MBF "caches", i.e. simple or structured arrays updated based on the current value of an associated
     * many body basis function object
     */
    struct Cache {
        /**
         * only used for debugging as a check that the current state of the cache corresponds to that of the MBF object
         * @return
         *  true if the current hash matches that of the last update
         */
        virtual bool is_valid() const = 0;
    protected:
        /**
         * should be updated only in the debug build any time an update is performed
         */
        defs::hash_t m_last_update_hash = 0;
    };

    /**
     * in most cases, a cache only requires a flat integer vector in which to store indices
     */
    struct SimpleContainer {
    protected:
        defs::inds m_inds;

    public:

        void clear();

        bool empty();
    };
};

#endif //M7_DECODEDMBF_H
