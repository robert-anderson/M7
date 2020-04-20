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
#include <src/defs.h>
#include "src/utils.h"

#ifdef HAVE_MPI
template<typename T>
static const MPI_Datatype mpi_type() {return MPI_Datatype();}

template<> const MPI_Datatype mpi_type<char>(){return MPI_CHAR;}
template<> const MPI_Datatype mpi_type<short int>(){return MPI_SHORT;}
template<> const MPI_Datatype mpi_type<int>(){return MPI_INT;}
template<> const MPI_Datatype mpi_type<long int>(){return MPI_LONG;}
template<> const MPI_Datatype mpi_type<long long int>(){return MPI_LONG_LONG_INT;}

template<> const MPI_Datatype mpi_type<unsigned char>(){return MPI_UNSIGNED_CHAR;}
template<> const MPI_Datatype mpi_type<unsigned short int>(){return MPI_UNSIGNED_SHORT;}
template<> const MPI_Datatype mpi_type<unsigned int>(){return MPI_UNSIGNED;}
template<> const MPI_Datatype mpi_type<unsigned long int>(){return MPI_UNSIGNED_LONG;}
template<> const MPI_Datatype mpi_type<unsigned long long int>(){return MPI_UNSIGNED_LONG_LONG;}

template<> const MPI_Datatype mpi_type<float>(){return MPI_FLOAT;}
template<> const MPI_Datatype mpi_type<double>(){return MPI_DOUBLE;}
template<> const MPI_Datatype mpi_type<long double>(){return MPI_LONG_DOUBLE;}

template<> const MPI_Datatype mpi_type<std::complex<float>>(){return MPI_FLOAT;}
template<> const MPI_Datatype mpi_type<std::complex<double>>(){return MPI_DOUBLE;}
template<> const MPI_Datatype mpi_type<std::complex<long double>>(){return MPI_LONG_DOUBLE;}

template<> const MPI_Datatype mpi_type<bool>(){return MPI_CXX_BOOL;}

const std::array<MPI_Op, 3> op_map {MPI_MAX, MPI_MIN, MPI_SUM};
#endif

struct mpi {

    static size_t nrank();

    static size_t irank();

    static std::string processor_name();

    static void barrier();

    enum MpiOp {MpiMax, MpiMin, MpiSum};

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
    static bool bcast(T &data) {
        return bcast(&data);
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

    static bool i_am(const size_t i);

    static bool i_am_root();

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
