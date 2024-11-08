//
// Created by Robert John Anderson on 08/02/2024.
//

#include "Owner.h"

Owner::Owner(uint_t i_rank_owner) : m_i(i_rank_owner){}

Owner Owner::local() {
    return {};
}

Owner Owner::shared(size_t i_rank) {
    return {i_rank};
}

Owner Owner::shared() {
    return shared(mpi::irank_world_shmem_root());
}

bool Owner::is_shared() const {
    return m_i < mpi::nrank();
}

bool Owner::is_local() const {
    return !(is_shared());
}

bool Owner::i_am_owner() const {
    return is_local() || m_i == mpi::irank();
}

uint_t Owner::irank_owner() const {
    return is_shared() ? m_i : mpi::irank();
}
