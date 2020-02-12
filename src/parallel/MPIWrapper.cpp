//
// Created by Robert John Anderson on 2020-02-07.
//

#include "MPIWrapper.h"

size_t mpi::nrank() {
#if USE_MPI
    int tmp;
    MPI_Comm_size(MPI_COMM_WORLD, &tmp);
    return tmp;
#else
    return 1;
#endif
}

size_t mpi::irank() {
#if USE_MPI
    int tmp;
    MPI_Comm_rank(MPI_COMM_WORLD, &tmp);
    return tmp;
#else
    return 0;
#endif
}

std::string mpi::processor_name() {
    int tmp;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &tmp);
    return std::string(processor_name, tmp);
}

void mpi::barrier() {
    MPI_Barrier(MPI_COMM_WORLD);
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
    if (!initialized()) MPI_Init(argc, argv);
}

void mpi::finalize() {
    if (initialized() && !finalized()) MPI_Finalize();
}

bool mpi::i_am(const size_t i) {
    return irank() == i;
}

bool mpi::i_am_root() {
    return i_am(0);
}
