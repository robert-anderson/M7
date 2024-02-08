//
// Created by Robert John Anderson on 08/02/2024.
//

#ifndef M7_OWNER_H
#define M7_OWNER_H

#include <M7_lib/defs.h>
#include <M7_lib/parallel/MPIWrapper.h>

/**
 * Relates to whether a Buffer's memory is shared or local. If shared, this class stores the world communicator rank
 * index that is designated the owner (i.e. the rank with the sole right to modify the buffer)
 */
class Owner {
    const uint_t m_i;
    // value of ~0ul means the buffer is to be private, not shared memory
    Owner(): m_i(~0ul){}
    Owner(uint_t i_rank_owner);
public:

    /**
     * Only this rank writes and reads the given buffer
     */
    static Owner local();

    /**
     * Only i_rank writes to the given buffer but all ranks in the same shared memory realm can read
     */
    static Owner shared(size_t i_rank);

    /**
     * Only the root rank in each shared memory realm writes to the given buffer but all ranks in the same shared
     * memory realm can read
     */
    static Owner shared();

    bool is_shared() const;

    bool is_local() const;

    bool i_am_owner() const;

    uint_t irank_owner() const;
};


#endif //M7_OWNER_H
