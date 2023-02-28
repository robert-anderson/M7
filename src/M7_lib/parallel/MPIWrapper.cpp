//
// Created by Robert John Anderson on 2020-02-07.
//


#include <utility>
#include "M7_lib/util/Integer.h"
#include "M7_lib/io/Logging.h"

#include "SharedArray.h"
#include "MPIWrapper.h"

void mpi::barrier() {
    MPI_Barrier(MPI_COMM_WORLD);
}

void mpi::barrier_on_node() {
    MPI_Barrier(g_node_comm);
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
    for (uint_t irank=0ul; irank<nrank(); ++irank) tmp.push_back(evenly_shared_count(nitem_global, irank));
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

bool mpi::i_am(uint_t i) {
    return irank() == i;
}

bool mpi::i_am_root() {
    return i_am(0);
}

bool mpi::on_node_i_am(uint_t i) {
    return irank_on_node() == i;
}

bool mpi::on_node_i_am_root() {
    return on_node_i_am(0);
}

bool mpi::is_node_root(uint_t irank) {
    return g_node_roots[irank];
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
    MPI_Comm_split_type(MPI_COMM_WORLD, OMPI_COMM_TYPE_CORE, irank(), MPI_INFO_NULL, &g_node_comm);
    MPI_Comm_size(g_node_comm, &tmp);
    g_nrank_on_node = tmp;
    MPI_Comm_rank(g_node_comm, &tmp);
    g_irank_on_node = tmp;
    mpi::all_gather(char(on_node_i_am_root()), g_node_roots);

    // only the node roots write their world communicator indices to the shared array
    g_my_node_root_irank = SharedScalar<uint_t>(mpi::irank());
}

void mpi::blocking_print(const str_t &str) {
    for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank) {
        if (mpi::i_am(irank)) {
            std::cout << str << std::endl;
        }
        mpi::barrier();
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

uint_t g_irank = 0;
uint_t g_nrank = 1;
str_t g_processor_name = "";
MPI_Comm g_node_comm;
uint_t g_irank_on_node = 0;
uint_t g_nrank_on_node = 1;
v_t<char> g_node_roots = {};
uint_t g_my_node_root_irank = 0;
int g_p2p_tag = 0;
