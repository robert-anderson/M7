//
// Created by Robert John Anderson on 2020-02-07.
//

#include "MPIWrapper.h"

#include <utility>
#include <M7_lib/io/Logging.h>

void mpi::barrier() {
    MPI_Barrier(MPI_COMM_WORLD);
}

void mpi::barrier_on_node() {
    MPI_Barrier(g_node_comm);
}

mpi::count_t mpi::evenly_shared_count(size_t nitem_global, size_t irank) {
    auto remainder = nitem_global % nrank();
    return nitem_global / nrank() + (irank < remainder);
}

mpi::count_t mpi::evenly_shared_count(size_t nitem_global) {
    return evenly_shared_count(nitem_global, mpi::irank());
}

mpi::counts_t mpi::evenly_shared_counts(size_t nitem_global) {
    mpi::counts_t tmp;
    tmp.reserve(nrank());
    for (size_t irank=0ul; irank<nrank(); ++irank) tmp.push_back(evenly_shared_count(nitem_global, irank));
    return tmp;
}

size_t mpi::evenly_shared_displ(size_t nitem_global, size_t irank) {
    auto ntb = std::min(irank, nitem_global % nrank());
    return ntb + irank * (nitem_global / nrank());
}

size_t mpi::evenly_shared_displ(size_t nitem_global) {
    return evenly_shared_displ(nitem_global, mpi::irank());
}

mpi::counts_t mpi::evenly_shared_displs(size_t nitem_global) {
    mpi::counts_t tmp;
    tmp.reserve(nrank());
    for (size_t irank=0ul; irank<nrank(); ++irank) tmp.push_back(evenly_shared_displ(nitem_global, irank));
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

bool mpi::i_am(const size_t& i) {
    return irank() == i;
}

bool mpi::i_am_root() {
    return i_am(0);
}

bool mpi::on_node_i_am(const size_t& i) {
    return irank_on_node() == i;
}

bool mpi::on_node_i_am_root() {
    return on_node_i_am(0);
}

void mpi::abort_(std::string message) {
    log::error_("Forcing MPI_Abort from this rank: {}", std::move(message));
    log::error_backtrace_();
    log::finalize();
    // SIGABRT is caught by IDEs for nice call stack debugging in the serial case
    if (mpi::nrank()==1) std::abort();
    MPI_Abort(MPI_COMM_WORLD, -1);
}

void mpi::abort(std::string message){
    if (mpi::nrank()==1)
        log::error("Reason: {}", std::move(message));
    else
        log::error_("Reason: {}", std::move(message));
    log::error_backtrace_();
    log::finalize();
    MPI_Barrier(MPI_COMM_WORLD);
    // SIGABRT is caught by IDEs for nice call stack debugging in the serial case
    if (mpi::nrank()==1) std::abort();
    MPI_Abort(MPI_COMM_WORLD, -1);
}


void mpi::setup_mpi_globals() {
    int tmp;
    MPI_Comm_size(MPI_COMM_WORLD, &tmp);
    g_nrank = tmp;
    ASSERT(g_nrank > 0)
    MPI_Comm_rank(MPI_COMM_WORLD, &tmp);
    g_irank = tmp;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &tmp);
    g_processor_name = std::string(processor_name, tmp);
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, irank(), MPI_INFO_NULL, &g_node_comm);
    MPI_Comm_size(g_node_comm, &tmp);
    g_nrank_on_node = tmp;
    MPI_Comm_rank(g_node_comm, &tmp);
    g_irank_on_node = tmp;
}

void mpi::blocking_print(const std::string &str) {
    for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) {
        if (mpi::i_am(irank)) {
            std::cout << str << std::endl;
        }
        mpi::barrier();
    }
}


size_t g_irank = 0;
size_t g_nrank = 1;
std::string g_processor_name = "";
MPI_Comm g_node_comm;
size_t g_irank_on_node = 0;
size_t g_nrank_on_node = 1;
int g_p2p_tag = 0;
