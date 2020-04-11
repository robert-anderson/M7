//
// Created by Robert John Anderson on 2020-02-07.
//

#include "MPIWrapper.h"

size_t mpi::nrank() {
#if HAVE_MPI
    int tmp;
    MPI_Comm_size(MPI_COMM_WORLD, &tmp);
    return tmp;
#else
    return 1;
#endif
}

size_t mpi::irank() {
#if HAVE_MPI
    int tmp;
    MPI_Comm_rank(MPI_COMM_WORLD, &tmp);
    return tmp;
#else
    return 0;
#endif
}

std::string mpi::processor_name() {
#if HAVE_MPI
    int tmp;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &tmp);
    return std::string(processor_name, tmp);
#endif
    return "node";
}

void mpi::barrier() {
#if HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
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
    if (!initialized()) MPI_Init(argc, argv);
#endif
}

void mpi::finalize() {
#if HAVE_MPI
    if (initialized() && !finalized()) MPI_Finalize();
#endif
}

bool mpi::i_am(const size_t i) {
    return irank() == i;
}

bool mpi::i_am_root() {
    return i_am(0);
}
