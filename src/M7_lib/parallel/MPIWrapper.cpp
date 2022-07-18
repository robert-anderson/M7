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

mpi::count_t mpi::evenly_shared_count(uint_t nitem_global, uint_t irank) {
    auto remainder = nitem_global % nrank();
    return nitem_global / nrank() + (irank < remainder);
}

mpi::count_t mpi::evenly_shared_count(uint_t nitem_global) {
    return evenly_shared_count(nitem_global, mpi::irank());
}

mpi::counts_t mpi::evenly_shared_counts(uint_t nitem_global) {
    mpi::counts_t tmp;
    tmp.reserve(nrank());
    for (uint_t irank=0ul; irank<nrank(); ++irank) tmp.push_back(evenly_shared_count(nitem_global, irank));
    return tmp;
}

uint_t mpi::evenly_shared_displ(uint_t nitem_global, uint_t irank) {
    auto ntb = std::min(irank, nitem_global % nrank());
    return ntb + irank * (nitem_global / nrank());
}

uint_t mpi::evenly_shared_displ(uint_t nitem_global) {
    return evenly_shared_displ(nitem_global, mpi::irank());
}

mpi::counts_t mpi::evenly_shared_displs(uint_t nitem_global) {
    mpi::counts_t tmp;
    tmp.reserve(nrank());
    for (uint_t irank=0ul; irank<nrank(); ++irank) tmp.push_back(evenly_shared_displ(nitem_global, irank));
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

bool mpi::i_am(const uint_t& i) {
    return irank() == i;
}

bool mpi::i_am_root() {
    return i_am(0);
}

bool mpi::on_node_i_am(const uint_t& i) {
    return irank_on_node() == i;
}

bool mpi::on_node_i_am_root() {
    return on_node_i_am(0);
}

void mpi::abort_(str_t message) {
    logging::error_("Forcing MPI_Abort from this rank: {}", std::move(message));
    //logging::error_backtrace_();
    logging::finalize();
    // SIGABRT is caught by IDEs for nice call stack debugging in the serial case
    if (mpi::nrank()==1) std::abort();
    MPI_Abort(MPI_COMM_WORLD, -1);
}

void mpi::abort(str_t message){
    if (mpi::nrank()==1)
        logging::error("Reason: {}", std::move(message));
    else
        logging::error_("Reason: {}", std::move(message));
    logging::error_backtrace_();
    logging::finalize();
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
    g_processor_name = str_t(processor_name, tmp);
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, irank(), MPI_INFO_NULL, &g_node_comm);
    MPI_Comm_size(g_node_comm, &tmp);
    g_nrank_on_node = tmp;
    MPI_Comm_rank(g_node_comm, &tmp);
    g_irank_on_node = tmp;
}

void mpi::blocking_print(const str_t &str) {
    for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank) {
        if (mpi::i_am(irank)) {
            std::cout << str << std::endl;
        }
        mpi::barrier();
    }
}


uint_t g_irank = 0;
uint_t g_nrank = 1;
str_t g_processor_name = "";
MPI_Comm g_node_comm;
uint_t g_irank_on_node = 0;
uint_t g_nrank_on_node = 1;
int g_p2p_tag = 0;
