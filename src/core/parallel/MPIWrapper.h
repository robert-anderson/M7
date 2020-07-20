//
// Created by Robert John Anderson on 2020-02-07.
//

#ifndef M7_MPIWRAPPER_H
#define M7_MPIWRAPPER_H

#include <array>

#ifdef HAVE_MPI

#include <mpi.h>

#endif

#include <iostream>
#include <cstring>
#include <src/core/util/defs.h>
#include "src/core/util/utils.h"

#ifdef HAVE_MPI

template<typename T>
static MPI_Datatype mpi_type() { return MPI_Datatype(); }

template<>
MPI_Datatype mpi_type<char>() { return MPI_CHAR; }

template<>
MPI_Datatype mpi_type<short int>() { return MPI_SHORT; }

template<>
MPI_Datatype mpi_type<int>() { return MPI_INT; }

template<>
MPI_Datatype mpi_type<long int>() { return MPI_LONG; }

template<>
MPI_Datatype mpi_type<long long int>() { return MPI_LONG_LONG_INT; }

template<>
MPI_Datatype mpi_type<unsigned char>() { return MPI_UNSIGNED_CHAR; }

template<>
MPI_Datatype mpi_type<unsigned short int>() { return MPI_UNSIGNED_SHORT; }

template<>
MPI_Datatype mpi_type<unsigned int>() { return MPI_UNSIGNED; }

template<>
MPI_Datatype mpi_type<unsigned long int>() { return MPI_UNSIGNED_LONG; }

template<>
MPI_Datatype mpi_type<unsigned long long int>() { return MPI_UNSIGNED_LONG_LONG; }

template<>
MPI_Datatype mpi_type<float>() { return MPI_FLOAT; }

template<>
MPI_Datatype mpi_type<double>() { return MPI_DOUBLE; }

template<>
MPI_Datatype mpi_type<long double>() { return MPI_LONG_DOUBLE; }

template<>
MPI_Datatype mpi_type<std::complex<float>>() { return MPI_COMPLEX; }

template<>
MPI_Datatype mpi_type<std::complex<double>>() { return MPI_DOUBLE_COMPLEX; }

template<>
MPI_Datatype mpi_type<std::complex<long double>>() { return MPI_CXX_LONG_DOUBLE_COMPLEX; }

template<>
MPI_Datatype mpi_type<bool>() { return MPI_CXX_BOOL; }


/*
 * For use with minloc and maxloc
 * MPI_FLOAT_INT: struct { float, int }
 * MPI_LONG_INT: struct { long, int }
 * MPI_DOUBLE_INT: struct { double, int }
 * MPI_SHORT_INT: struct { short, int }
 * MPI_2INT: struct { int, int }
 * MPI_LONG_DOUBLE_INT: struct { long double, int }
 */

template<typename T>
static MPI_Datatype mpi_pair_type() { return MPI_Datatype(); }

template<>
MPI_Datatype mpi_pair_type<float>() {return MPI_FLOAT_INT;}

template<>
MPI_Datatype mpi_pair_type<long>() {return MPI_LONG_INT;}

template<>
MPI_Datatype mpi_pair_type<double>() {return MPI_DOUBLE_INT;}

template<>
MPI_Datatype mpi_pair_type<short>() {return MPI_SHORT_INT;}

template<>
MPI_Datatype mpi_pair_type<int>() {return MPI_2INT;}

template<>
MPI_Datatype mpi_pair_type<long double>() {return MPI_LONG_DOUBLE_INT;}


const std::array<MPI_Op, 5> op_map{MPI_MAX, MPI_MIN, MPI_SUM, MPI_LAND, MPI_LOR};
const std::array<MPI_Op, 2> pair_op_map{MPI_MAXLOC, MPI_MINLOC};

enum MpiOp {
    MpiMax, MpiMin, MpiSum, MpiLand, MpiLor
};

enum MpiPairOp {
    MpiMaxLoc, MpiMinLoc
};


#endif

extern size_t g_irank;
extern size_t g_nrank;
extern std::string g_processor_name;
#ifdef HAVE_MPI
extern MPI_Comm g_node_comm;
#endif
extern size_t g_irank_on_node;
extern size_t g_nrank_on_node;

struct mpi {
    static void setup_mpi_globals(){
#ifdef HAVE_MPI
        int tmp;
        MPI_Comm_size(MPI_COMM_WORLD, &tmp);
        g_nrank = tmp;
        ASSERT(g_nrank>0)
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

    static const size_t &nrank(){
        return g_nrank;
    }

    static const size_t &irank(){
        return g_irank;
    }

    static const size_t &nrank_on_node(){
        return g_nrank_on_node;
    }

    static const size_t &irank_on_node(){
        return g_irank_on_node;
    }

    static const std::string &processor_name(){
        return g_processor_name;
    }

    static MPI_Comm* node_communicator(){
#ifdef HAVE_MPI
        return &g_node_comm;
#else
        return nullptr;
#endif
    }

    static void barrier();
    static void barrier_on_node();

/*
 * int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm)
 */

private:
    template<typename T>
    static bool reduce(const T *send, T *recv, MpiOp op, size_t ndata = 1, size_t iroot = 0) {
#ifdef HAVE_MPI
        return MPI_Reduce(send, recv, ndata, mpi_type<T>(), op_map[op], iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        std::memcpy(recv, send, sizeof(T)*ndata);
        return true;
#endif
    }

    template<typename T>
    static bool all_reduce(const T *send, T *recv, MpiOp op, size_t ndata = 1) {
#ifdef HAVE_MPI
        return MPI_Allreduce(send, recv, ndata, mpi_type<T>(), op_map[op], MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        std::memcpy(recv, send, sizeof(T)*ndata);
        return true;
#endif
    }

    template<typename T>
    static bool all_reduce(const std::pair<T, int> *send, std::pair<T, int> *recv, MpiPairOp op, size_t ndata = 1) {
#ifdef HAVE_MPI
        return MPI_Allreduce(send, recv, ndata, mpi_pair_type<T>(), pair_op_map[op], MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        std::memcpy(recv, send, sizeof(T)*ndata);
        return true;
#endif
    }

public:
    /*
     * MAX REDUCE CONVENIENCE METHODS
     */
    template<typename T>
    static bool max(const T *send, T *recv, size_t ndata = 1, size_t iroot = 0) {
        return reduce(send, recv, MpiMax, ndata, iroot);
    }

    template<typename T>
    static bool all_max(const T *send, T *recv, size_t ndata = 1) {
        return all_reduce(send, recv, MpiMax, ndata);
    }

    template<typename T>
    static T max(const T *send, size_t iroot = 0) {
        T recv;
        max(send, &recv, 1, iroot);
        return recv;
    }

    template<typename T>
    static T all_max(const T *send) {
        T recv;
        all_max(send, &recv, 1);
        return recv;
    }

    template<typename T>
    static T all_max(const T &send) {
        return all_max(&send);
    }

    /*
     * MIN REDUCE CONVENIENCE METHODS
     */
    template<typename T>
    static bool min(const T *send, T *recv, size_t ndata = 1, size_t iroot = 0) {
        return reduce(send, recv, MpiMin, ndata, iroot);
    }

    template<typename T>
    static bool all_min(const T *send, T *recv, size_t ndata = 1) {
        return all_reduce(send, recv, MpiMin, ndata);
    }

    template<typename T>
    static T min(const T *send, size_t iroot = 0) {
        T recv;
        min(send, &recv, 1, iroot);
        return recv;
    }

    template<typename T>
    static T all_min(const T *send) {
        T recv;
        all_min(send, &recv, 1);
        return recv;
    }

    template<typename T>
    static T all_min(const T &send) {
        return all_min(&send);
    }

    /*
     * SUM REDUCE CONVENIENCE METHODS
     */
    template<typename T>
    static bool sum(const T *send, T *recv, size_t ndata = 1, size_t iroot = 0) {
        return reduce(send, recv, MpiSum, ndata, iroot);
    }

    template<typename T>
    static bool all_sum(const T *send, T *recv, size_t ndata = 1) {
        return all_reduce(send, recv, MpiSum, ndata);
    }

    template<typename T>
    static T sum(const T *send, size_t iroot = 0) {
        T recv;
        sum(send, &recv, 1, iroot);
        return recv;
    }

    template<typename T>
    static T all_sum(const T *send) {
        T recv;
        all_sum(send, &recv, 1);
        return recv;
    }

    template<typename T>
    static T all_sum(const T &send) {
        return all_sum(&send);
    }

    /*
     * LOGICAL-AND REDUCE CONVENIENCE FUNCTIONS
     */

    template<typename T>
    static bool land(const T *send, T *recv, size_t ndata = 1, size_t iroot = 0) {
        return reduce(send, recv, MpiLand, ndata, iroot);
    }

    template<typename T>
    static bool all_land(const T *send, T *recv, size_t ndata = 1) {
        return all_reduce(send, recv, MpiLand, ndata);
    }

    template<typename T>
    static T land(const T *send, size_t iroot = 0) {
        T recv;
        land(send, &recv, 1, iroot);
        return recv;
    }

    template<typename T>
    static T all_land(const T *send) {
        T recv;
        all_land(send, &recv, 1);
        return recv;
    }

    template<typename T>
    static T all_land(const T &send) {
        return all_land(&send);
    }


    /*
     * BCAST
     */
    template<typename T>
    static bool bcast(T *data, const size_t ndata = 1, size_t iroot = 0) {
#ifdef HAVE_MPI
        return MPI_Bcast((void *) data, ndata, mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        return true;
#endif
    }

    template<typename T>
    static bool bcast(T &data, size_t iroot = 0) {
        return bcast(&data, 1, iroot);
    }

    template<typename T>
    static bool bcast(std::vector<T> &data, size_t ndata = 0, size_t iroot = 0) {
#ifdef HAVE_MPI
        if (!ndata) ndata = data.size();
        return MPI_Bcast((void *) data.data(), ndata, mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        return true;
#endif
    }


    /*
     * MAXLOC CONVENIENCE FUNCTIONS
     */

    template<typename T>
    static bool all_maxloc(const T &send, std::pair<T, int> &recv) {
        std::pair<T, int> tmp {send, irank()};
        return all_reduce(&tmp, &recv, MpiMaxLoc, 1);
    }

    /*
     * MINLOC CONVENIENCE FUNCTIONS
     */

    template<typename T>
    static bool all_minloc(const T &send, std::pair<T, int> &recv) {
        std::pair<T, int> tmp {send, irank()};
        return all_reduce(&tmp, &recv, MpiMinLoc, 1);
    }


    template<typename T>
    static bool all_to_all(const T *send, const size_t nsend, T *recv, const size_t nrecv) {
#ifdef HAVE_MPI
        return MPI_Alltoall((void *) send, nsend, mpi_type<T>(),
                            (void *) recv, nrecv, mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        ASSERT(nsend == nrecv);
        return all_reduce(send, recv, MpiMax, nsend);
#endif
    }

    template<typename T>
    static bool all_to_all(const std::vector<T> &send, std::vector<T> &recv) {
        return all_to_all(send.data(), 1, recv.data(), 1);
    }

private:
    /*
     * we don't expose this interfacing method, because defs::inds are always std::vector<size_t>
     */
    template<typename T>
    static bool all_to_allv(
            const T *send, const int *sendcounts, const int *senddispls,
            T *recv, const int *recvcounts, const int *recvdispls) {
#ifdef HAVE_MPI
        return MPI_Alltoallv(
                (void *) send, sendcounts, senddispls, mpi_type<T>(),
                (void *) recv, recvcounts, recvdispls, mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        return all_to_all(send, senddispls[0]+sendcounts[0], recv, recvdispls[0]+recvcounts[0]);
#endif
    }

    template<typename T>
    static bool all_gatherv(
            const T *send, const int sendcount,
            T *recv, const int *recvcounts, const int *displs) {
#ifdef HAVE_MPI
        return MPI_Allgatherv(
                (void *) send, sendcount, mpi_type<T>(),
                (void *) recv, recvcounts, displs, mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        return all_to_all(send, sendcount, recv, recvcounts[0]);
#endif
    }

public:
    template<typename T>
    static bool all_to_allv(
            const T *send, const defs::inds &sendcounts, const defs::inds &senddispls,
            T *recv, const defs::inds &recvcounts, const defs::inds &recvdispls) {
        return all_to_allv(send, std::vector<int>(sendcounts.begin(), sendcounts.end()).data(),
                           std::vector<int>(senddispls.begin(), senddispls.end()).data(),
                           recv, std::vector<int>(recvcounts.begin(), recvcounts.end()).data(),
                           std::vector<int>(recvdispls.begin(), recvdispls.end()).data());
    }

    template<typename T>
    static bool all_gatherv(
            const T *send, const size_t &sendcount,
            T *recv, const defs::inds &recvcounts, const defs::inds &recvdispls) {
        return all_gatherv(send, sendcount,
                           recv, std::vector<int>(recvcounts.begin(), recvcounts.end()).data(),
                           std::vector<int>(recvdispls.begin(), recvdispls.end()).data());
    }

    template<typename T>
    static bool all_gather(
            const T *send, const int sendcount, T *recv, const int recvcount) {
#ifdef HAVE_MPI
        return MPI_Allgather(
                (void *) send, sendcount, mpi_type<T>(), recv, recvcount, mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        return all_to_all(send, sendcount, recv, recvcount);
#endif
    }

    static bool i_am(const size_t i);

    static bool on_node_i_am(const size_t i);

    static bool i_am_root();

    static bool on_node_i_am_root();

    static void rank_print(const std::string s, size_t irank) {
        if (i_am(irank)) std::cout << "rank " << irank << ": " << s << std::endl;
        barrier();
    }

    static void root_print(const std::string s) {
        rank_print(s, 0);
    }

    static bool initialized();

    static bool finalized();

    static void initialize(int *argc = NULL, char ***argv = NULL);

    static void finalize();

};

#endif //M7_MPIWRAPPER_H
