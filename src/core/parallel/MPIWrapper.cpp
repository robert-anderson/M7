//
// Created by Robert John Anderson on 2020-02-07.
//

#include "MPIWrapper.h"
#include "src/core/io/Logging.h"

void mpi::barrier() {
#if HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void mpi::barrier_on_node() {
#if HAVE_MPI
    MPI_Barrier(g_node_comm);
#endif
}

bool mpi::initialized() {
#if HAVE_MPI
    int tmp;
    MPI_Initialized(&tmp);
    return tmp;
#endif
    return true;
}

bool mpi::finalized() {
#if HAVE_MPI
    int tmp;
    MPI_Finalized(&tmp);
    return tmp;
#endif
    return true;
}

void mpi::initialize(int *argc, char ***argv) {
#if HAVE_MPI
    if (!initialized()) {
        MPI_Init(argc, argv);
        setup_mpi_globals();
    }
#endif
}

void mpi::finalize() {
#if HAVE_MPI
    if (initialized() && !finalized()) {
        MPI_Finalize();
    }
#endif
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
#ifdef HAVE_MPI
    log::error_("Aborting all MPI processes: \"{}\"", message);
    log::finalize();
    // SIGABRT is caught by IDEs for nice call stack debugging in the serial case
    if (mpi::nrank()==1) std::abort();
    MPI_Abort(MPI_COMM_WORLD, -1);
#else
    std::cout << "Aborting: \"" << message << "\""<< std::endl;
        exit(0);
#endif
}

void mpi::abort(std::string message, float f){
#ifdef HAVE_MPI
    log::error_("Aborting all MPI processes: \"{}\"", message);
    log::finalize();
    MPI_Barrier(MPI_COMM_WORLD);
    // SIGABRT is caught by IDEs for nice call stack debugging in the serial case
    if (mpi::nrank()==1) std::abort();
    MPI_Abort(MPI_COMM_WORLD, -1);
#else
    std::cout << "Aborting: \"" << message << "\""<< std::endl;
        exit(0);
#endif
}


void mpi::setup_mpi_globals() {
#ifdef HAVE_MPI
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
#endif
}


size_t g_irank = 0;
size_t g_nrank = 1;
std::string g_processor_name = "";
#ifdef HAVE_MPI
MPI_Comm g_node_comm;
#endif
size_t g_irank_on_node = 0;
size_t g_nrank_on_node = 1;