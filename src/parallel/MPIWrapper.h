//
// Created by Robert John Anderson on 2020-02-07.
//

#ifndef M7_MPIWRAPPER_H
#define M7_MPIWRAPPER_H

#include <array>
#include <mpi.h>
#include "../data/DataTable.h"

#define USE_MPI 1

static const std::array<MPI_Datatype, nnumeric> type_to_mpi_type = {
        (MPI_Datatype) MPI_COMPLEX,
        (MPI_Datatype) MPI_DOUBLE_COMPLEX,
        (MPI_Datatype) MPI_CXX_LONG_DOUBLE_COMPLEX,
        (MPI_Datatype) MPI_FLOAT,
        (MPI_Datatype) MPI_DOUBLE,
        (MPI_Datatype) MPI_LONG_DOUBLE,
        (MPI_Datatype) MPI_CHAR,
        (MPI_Datatype) MPI_SHORT,
        (MPI_Datatype) MPI_INT,
        (MPI_Datatype) MPI_LONG_INT,
        (MPI_Datatype) MPI_LONG_LONG_INT,
        (MPI_Datatype) MPI_UNSIGNED_CHAR,
        (MPI_Datatype) MPI_UNSIGNED_SHORT,
        (MPI_Datatype) MPI_UNSIGNED,
        (MPI_Datatype) MPI_UNSIGNED_LONG,
        (MPI_Datatype) MPI_UNSIGNED_LONG_LONG,
        (MPI_Datatype) MPI_CXX_BOOL
};

template<typename T>
MPI_Datatype mpi_type() {
    return type_to_mpi_type[type_number<T>];
}

class MPIWrapper {
    size_t m_nrank = 1;
    size_t m_irank = 0;
    std::string m_processor_name = "";
public:
    MPIWrapper() {
#if USE_MPI
        int initialized;
        MPI_Initialized(&initialized);
        assert(initialized);
        int tmp;
        MPI_Comm_size(MPI_COMM_WORLD, &tmp);
        m_nrank = tmp;
        MPI_Comm_rank(MPI_COMM_WORLD, &tmp);
        m_irank = tmp;
        char processor_name[MPI_MAX_PROCESSOR_NAME];
        MPI_Get_processor_name(processor_name, &tmp);
        m_processor_name = std::string(processor_name, tmp);
#endif
    }

    void barrier() {
        MPI_Barrier(MPI_COMM_WORLD);
    }
/*
 * int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm)
 */
private:
    template<typename T>
    bool reduce(const T *send, T *recv, MPI_Op op, size_t ndata = 1, size_t iroot = 0) {
        return MPI_Reduce(send, recv, ndata, mpi_type<T>(), op, iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    bool all_reduce(const T *send, T *recv, MPI_Op op, size_t ndata = 1) {
        return MPI_Allreduce(send, recv, ndata, mpi_type<T>(), op, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

public:

    /*
     * MAX REDUCE CONVENIENCE METHODS
     */
    template<typename T>
    bool max(const T *send, T *recv, size_t ndata = 1, size_t iroot = 0) {
        return reduce(send, recv, MPI_MAX, ndata, iroot);
    }

    template<typename T>
    bool all_max(const T *send, T *recv, size_t ndata = 1) {
        return all_reduce(send, recv, MPI_MAX, ndata);
    }

    template<typename T>
    T max(const T *send, size_t iroot = 0) {
        T *recv;
        max(send, recv, 1, iroot);
        return *recv;
    }

    template<typename T>
    T all_max(const T *send) {
        T *recv;
        all_max(send, recv, 1);
        return *recv;
    }

    template<typename T>
    T all_max(const T &send) {
        return all_max(&send);
    }

    /*
     * MIN REDUCE CONVENIENCE METHODS
     */
    template<typename T>
    bool min(const T *send, T *recv, size_t ndata = 1, size_t iroot = 0) {
        return reduce(send, recv, MPI_MIN, ndata, iroot);
    }

    template<typename T>
    bool all_min(const T *send, T *recv, size_t ndata = 1) {
        return all_reduce(send, recv, MPI_MIN, ndata);
    }

    template<typename T>
    T min(const T *send, size_t iroot = 0) {
        T *recv;
        min(send, recv, 1, iroot);
        return *recv;
    }

    template<typename T>
    T all_min(const T *send) {
        T *recv;
        all_min(send, recv, 1);
        return *recv;
    }

    template<typename T>
    T all_min(const T &send) {
        return all_min(&send);
    }

    /*
     * SUM REDUCE CONVENIENCE METHODS
     */
    template<typename T>
    bool sum(const T *send, T *recv, size_t ndata = 1, size_t iroot = 0) {
        return reduce(send, recv, MPI_SUM, ndata, iroot);
    }

    template<typename T>
    bool all_sum(const T *send, T *recv, size_t ndata = 1) {
        return all_reduce(send, recv, MPI_SUM, ndata);
    }

    template<typename T>
    T sum(const T *send, size_t iroot = 0) {
        T *recv;
        sum(send, recv, 1, iroot);
        return *recv;
    }

    template<typename T>
    T all_sum(const T *send) {
        T *recv;
        all_sum(send, recv, 1);
        return *recv;
    }

    template<typename T>
    T all_sum(const T &send) {
        return all_sum(&send);
    }

    /*
     * BCAST
     */
    template<typename T>
    bool bcast(T *data, const size_t ndata = 1, size_t iroot = 0) {
        return MPI_Bcast((void *) data, ndata, mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    bool bcast(T &data) {
        return bcast(&data);
    }

    template<typename T>
    bool bcast(std::vector<T> &data, size_t ndata = 0, size_t iroot = 0) {
        if (!ndata) ndata = data.size();
        return MPI_Bcast((void *) data.data(), ndata, mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    bool all_to_all(const T *send, const size_t nsend, T *recv, const size_t nrecv) {
        return MPI_Alltoall((void *) send, nsend, mpi_type<T>(),
                            (void *) recv, nrecv, mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    bool all_to_all(const std::vector<T> &send, std::vector<T> &recv) {
        return all_to_all(send.data(), 1, recv.data(), 1);
    }

private:
    /*
     * we don't expose this interfacing method, because defs::inds are always std::vector<size_t>
     */
    template<typename T>
    bool all_to_allv(
            const T *send, const int *sendcounts, const int *senddispls,
            T *recv, const int *recvcounts, const int *recvdispls) {
        return MPI_Alltoallv(
                (void *) send, sendcounts, senddispls, mpi_type<T>(),
                (void *) recv, recvcounts, recvdispls, mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
    }

public:

    template<typename T>
    bool all_to_allv(
            const T *send, const defs::inds &sendcounts, const defs::inds &senddispls,
            T *recv, const defs::inds &recvcounts, const defs::inds &recvdispls) {
        return all_to_allv(send, std::vector<int>(sendcounts.begin(), sendcounts.end()).data(),
                           std::vector<int>(senddispls.begin(), senddispls.end()).data(),
                           recv, std::vector<int>(recvcounts.begin(), recvcounts.end()).data(),
                           std::vector<int>(recvdispls.begin(), recvdispls.end()).data());
    }

    size_t nrank() const {
        return m_nrank;
    }

    size_t irank() const {
        return m_irank;
    }

    const std::string &processor_name() const {
        return m_processor_name;
    }

    bool i_am(const size_t i) {
        return irank() == i;
    }

    bool i_am_root() {
        return i_am(0);
    }

    void rank_print(const std::string s, size_t irank) {
        if (i_am(irank)) std::cout << "rank " << irank << ": " << s << std::endl;
        barrier();
    }

    void root_print(const std::string s) {
        rank_print(s, 0);
    }

    static void initialize(int *argc = NULL, char ***argv = NULL) {
        int initialized;
        MPI_Initialized(&initialized);
        if (!initialized) MPI_Init(argc, argv);
    }

    static void finalize() {
        int finalized;
        MPI_Finalized(&finalized);
        if (!finalized) MPI_Finalize();
    }

};

#endif //M7_MPIWRAPPER_H
