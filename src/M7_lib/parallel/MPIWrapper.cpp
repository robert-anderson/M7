//
// Created by Robert John Anderson on 2020-02-07.
//


#include <utility>
#include <numeric>
#include <set>
#include "M7_lib/util/Integer.h"
#include "M7_lib/io/Logging.h"

#include "SharedArray.h"
#include "MPIWrapper.h"

void mpi::barrier(mpi::Realm realm) {
    realm == World ? MPI_Barrier(MPI_COMM_WORLD) : MPI_Barrier(g_shmem_comm);
}

mpi::count_t mpi::evenly_shared_count(uint_t nitem_global, uint_t irank) {
    return integer::evenly_shared_count(nitem_global, irank, nrank());
}

mpi::count_t mpi::evenly_shared_count(uint_t nitem_global) {
    return evenly_shared_count(nitem_global, mpi::irank());
}

mpi::countv_t mpi::evenly_shared_counts(uint_t nitem_global) {
    mpi::countv_t tmp;
    tmp.reserve(nrank());
    for (uint_t irank = 0ul; irank < nrank(); ++irank) tmp.push_back(evenly_shared_count(nitem_global, irank));
    return tmp;
}

uint_t mpi::evenly_shared_displ(uint_t nitem_global, uint_t irank) {
    return integer::evenly_shared_offset(nitem_global, irank, nrank());
}

uint_t mpi::evenly_shared_displ(uint_t nitem_global) {
    return evenly_shared_displ(nitem_global, mpi::irank());
}

mpi::countv_t mpi::evenly_shared_displs(uint_t nitem_global) {
    mpi::countv_t tmp;
    tmp.reserve(nrank());
    for (uint_t irank = 0ul; irank < nrank(); ++irank) tmp.push_back(evenly_shared_displ(nitem_global, irank));
    return tmp;
}

bool mpi::initialized() {
    int tmp;
    MPI_Initialized(&tmp);
    return tmp;
}

bool mpi::finalized() {
    int tmp;
    MPI_Finalized(&tmp);
    return tmp;
}

void mpi::initialize(int *argc, char ***argv) {
    if (!initialized()) {
        MPI_Init(argc, argv);
        setup_mpi_globals();
    }
}

void mpi::finalize() {
    if (initialized() && !finalized()) {
        MPI_Finalize();
    }
}

bool mpi::i_am(uint_t i, Realm realm) {
    return irank(realm) == i;
}

bool mpi::i_am_root(Realm realm) {
    return i_am(0, realm);
}

uint_t mpi::nshmem() {
    return g_iranks_world_in_shmem_realms.size();
}

bool mpi::is_root(uint_t irank_world, Realm realm) {
    return realm == World ? !irank_world : !g_iranks_shmem[irank_world];
}

void mpi::abort_(str_t message) {
    logging::error_("Forcing MPI_Abort from this rank: {}", std::move(message));
    //logging::error_backtrace_();
    logging::finalize();
    // SIGABRT is caught by IDEs for nice call stack debugging in the serial case
    if (mpi::nrank() == 1) std::abort();
    MPI_Abort(MPI_COMM_WORLD, -1);
}

void mpi::abort(str_t message){
    if (mpi::nrank() == 1)
        logging::error("Reason: {}", std::move(message));
    else
        logging::error_("Reason: {}", std::move(message));
    logging::error_backtrace_();
    logging::finalize();
    MPI_Barrier(MPI_COMM_WORLD);
    // SIGABRT is caught by IDEs for nice call stack debugging in the serial case
    if (mpi::nrank() == 1) std::abort();
    MPI_Abort(MPI_COMM_WORLD, -1);
}


void mpi::setup_mpi_globals() {
    int tmp;
    g_world_comm = MPI_COMM_WORLD;
    // get the size of the world communicator (i.e. the total number of MPI ranks)
    MPI_Comm_size(g_world_comm, &tmp);
    g_nrank_world = tmp;
    ASSERT(tmp > 0);
    // get the index of this rank in the world communicator
    MPI_Comm_rank(g_world_comm, &tmp);
    g_irank_world = tmp;
    // split the world communicator by shared memory region into sub-communicators
    MPI_Comm_split_type(g_world_comm, MPI_COMM_TYPE_SHARED, irank(), MPI_INFO_NULL, &g_shmem_comm);
    // get the size of the shared memory communicator in which this rank resides
    MPI_Comm_size(g_shmem_comm, &tmp);
    g_nrank_shmem = tmp;
    // get the index of this rank in the shared memory communicator in which it resides
    MPI_Comm_rank(g_shmem_comm, &tmp);
    g_irank_shmem = tmp;

    {
        // make a group containing all ranks in the world communicator
        MPI_Group world_group;
        MPI_Comm_group(MPI_COMM_WORLD, &world_group);
        // make a group containing all ranks in the shared memory communicator
        MPI_Group shmem_group;
        MPI_Comm_group(g_shmem_comm, &shmem_group);
        std::vector<int> world_iranks(nrank());
        std::iota(world_iranks.begin(), world_iranks.end(), 0);
        std::vector<int> shmem_iranks(nrank());
        // translate the world rank indices into shared memory rank indices
        MPI_Group_translate_ranks(world_group, nrank(), world_iranks.data(), shmem_group, shmem_iranks.data());
        // find the world index of the root of the shared memory communicator in which this rank resides
        auto it = std::find(shmem_iranks.cbegin(), shmem_iranks.cend(), 0);
        // the root should have been found somewhere
        assert(it != shmem_iranks.cend());
        const uint_t shmem_root_irank = std::distance(shmem_iranks.cbegin(), it);
        // gather all into the global vector mapping all ranks to their root rank in shared memory
        mpi::all_gather(shmem_root_irank, g_shmem_root_iranks_world);
    }

    {
        // use map to get a unique, ordered vector of realms with their ranks
        std::map<uint_t, std::set<uint_t>> tmp_map;
        for (uint_t irank = 0ul; irank < nrank(); ++irank) {
            const auto irank_shmem_root = g_shmem_root_iranks_world[irank];
            auto it = tmp_map.find(irank_shmem_root);
            if (it == tmp_map.end()) it = tmp_map.insert({irank_shmem_root, {}}).first;
            it->second.insert(irank);
        }
        g_ishmems.resize(nrank());
        uint_t ishmem = 0ul;
        for (auto pair: tmp_map) {
            g_irank_root_in_shmem_realms.push_back(pair.first);
            g_iranks_world_in_shmem_realms.emplace_back(pair.second.cbegin(), pair.second.cend());
            for (auto irank: g_iranks_world_in_shmem_realms.back()) g_ishmems[irank] = ishmem;
            g_nrank_in_shmem_realms.push_back(g_iranks_world_in_shmem_realms.back().size());
            ++ishmem;
        }
    }
}

void mpi::filter(bool cond, uintv_t& iranks) {
    /*
     * use v_t<unsigned char> as container, since the std::vector<bool> can't be used due to bitfield optimisation
     */
    v_t<unsigned char> gathered;
    const unsigned char c = cond;
    mpi::all_gather(c, gathered);
    iranks.clear();
    for (uint_t i = 0ul; i < gathered.size(); ++i) if (gathered[i]) iranks.push_back(i);
}

uintv_t mpi::filter(bool cond) {
    uintv_t ranks;
    filter(cond, ranks);
    return ranks;
}

MPI_Comm mpi::g_world_comm;
uint_t mpi::g_irank_world = 0;
uint_t mpi::g_nrank_world = 1;
MPI_Comm mpi::g_shmem_comm;
uint_t mpi::g_irank_shmem = 0;
uint_t mpi::g_nrank_shmem = 1;
uintv_t mpi::g_iranks_shmem = {};
uintv_t mpi::g_shmem_root_iranks_world = {};
uintv_t mpi::g_irank_root_in_shmem_realms = {};
v_t<uintv_t> mpi::g_iranks_world_in_shmem_realms = {};
uintv_t mpi::g_nrank_in_shmem_realms = {};
uintv_t mpi::g_ishmems = {};

int mpi::g_p2p_tag = 0;
std::array<v_t<buf_t>, mpi::g_types.size()> mpi::g_send_reduce_buffers = {};
std::array<v_t<buf_t>, mpi::g_types.size()> mpi::g_recv_reduce_buffers = {};
