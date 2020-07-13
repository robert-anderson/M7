//
// Created by Robert John Anderson on 2020-02-07.
//

#include "MPIWrapper.h"


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
    if (!initialized()) {
        MPI_Init(argc, argv);
        setup_mpi_globals();
    }
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

bool mpi::on_node_i_am(const size_t i) {
    return irank_on_node() == i;
}

bool mpi::on_node_i_am_root() {
    return on_node_i_am(0);
}
