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

/**
 * With MPI, we potentially have a typing issue. By default MPI libraries are typically
 * compiled with the displacement and count types set to 32-bit signed ints. This is
 * problematic whenever we want to communicate blocks of more than 2147483647 (2^32-1)
 * words (over 17GB worth of doubles). M7 uses unsigned longs (size_t) throughout for
 * integer offsets, indices, and sizes - a mismatch this wrapper must gracefully resolve.
 *
 * The utils::safe_narrow methods will raise a runtime exception if a narrowing conversion
 * would result in a loss of information if the SAFE_NARROWING macro is defined. In normal
 * usage, such an overflow is fairly unlikely. But we do not leave it to chance that such a
 * large number of words will never be communicated in practice.
 *
 * defs::mpi_count and defs::mpi_counts are the scalar and vector typedefs for the MPI integer
 *
 * The cleanest way to acheive this aim is to expose only the overloads that use size_t and
 * std::vector<size_t> (defs::inds) to define counts and displs. if sizeof(defs::mpi_count)
 * is less than sizeof(size_t), a narrowing conversion is required, but since this is at the
 * level of communication, and not inside a main loop, this minor overhead will more than
 * pay for itself in code clarity. Hence, **defs::mpi_count and defs::mpi_counts should never
 * be seen outside this module**.
 *
 * It is worth pointing out that the MPI shared memory windows employed in SharedArray.h and
 * its dependents use the MPI_Aint type to store sizes in chars. This size is large enough
 * to address any byte of memory, and so worries about narrowing in the MPI interface may be
 * forgotten in this case.
 */

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
MPI_Datatype mpi_pair_type<float>() { return MPI_FLOAT_INT; }

template<>
MPI_Datatype mpi_pair_type<long>() { return MPI_LONG_INT; }

template<>
MPI_Datatype mpi_pair_type<double>() { return MPI_DOUBLE_INT; }

template<>
MPI_Datatype mpi_pair_type<short>() { return MPI_SHORT_INT; }

template<>
MPI_Datatype mpi_pair_type<int>() { return MPI_2INT; }

template<>
MPI_Datatype mpi_pair_type<long double>() { return MPI_LONG_DOUBLE_INT; }


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

    static void setup_mpi_globals() {
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

    static const size_t &nrank() {
        return g_nrank;
    }

    static const size_t &irank() {
        return g_irank;
    }

    static const size_t &nrank_on_node() {
        return g_nrank_on_node;
    }

    static const size_t &irank_on_node() {
        return g_irank_on_node;
    }

    static const std::string &processor_name() {
        return g_processor_name;
    }

    static MPI_Comm *node_communicator() {
#ifdef HAVE_MPI
        return &g_node_comm;
#else
        return nullptr;
#endif
    }

    static void barrier();

    static void barrier_on_node();

    static defs::mpi_count snrw(const size_t& i){
        return utils::safe_narrow<defs::mpi_count>(i);
    }

    static defs::mpi_counts snrw(const std::vector<size_t>& v){
        return utils::safe_narrow<defs::mpi_count>(v);
    }

    template<typename T>
    static void counts_to_displs_consec(const std::vector<T> &sizes, std::vector<T> &displs) {
        ASSERT(sizes.size() == displs.size())
        displs[0] = 0;
        for (size_t i = 1ul; i < sizes.size(); ++i) displs[i] = displs[i - 1] + sizes[i - 1];
    }


private:
    template<typename T>
    static bool reduce(const T *send, T *recv, MpiOp op, size_t ndata = 1, size_t iroot = 0) {
#ifdef HAVE_MPI
        return MPI_Reduce(send, recv, snrw(ndata), mpi_type<T>(), op_map[op], snrw(iroot), MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        std::memcpy(recv, send, sizeof(T)*ndata);
        return true;
#endif
    }

    template<typename T>
    static bool all_reduce(const T *send, T *recv, MpiOp op, size_t ndata = 1) {
#ifdef HAVE_MPI
        return MPI_Allreduce(send, recv, snrw(ndata), mpi_type<T>(), op_map[op], MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        std::memcpy(recv, send, sizeof(T)*ndata);
        return true;
#endif
    }

    template<typename T>
    static bool all_reduce(const std::pair<T, size_t> *send, std::pair<T, size_t> *recv,
                           MpiPairOp op, size_t ndata = 1) {
#ifdef HAVE_MPI
        const std::pair<T, defs::mpi_count> tmp_send{send->first, snrw(send->second)};
        std::pair<T, defs::mpi_count> tmp_recv;
        auto res = MPI_Allreduce(&tmp_send, &tmp_recv, ndata, mpi_pair_type<T>(), pair_op_map[op], MPI_COMM_WORLD) == MPI_SUCCESS;
        recv->first = tmp_recv.first;
        recv->second = tmp_recv.second;
        return res;
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
    static bool bcast(T *data, size_t ndata = 1, size_t iroot = 0) {
#ifdef HAVE_MPI
        return MPI_Bcast((void *) data, snrw(ndata), mpi_type<T>(), snrw(iroot), MPI_COMM_WORLD) == MPI_SUCCESS;
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
        return MPI_Bcast((void *) data.data(), snrw(ndata), mpi_type<T>(), snrw(iroot), MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        return true;
#endif
    }


    /*
     * MAXLOC CONVENIENCE FUNCTIONS
     */

    template<typename T>
    static bool all_maxloc(const T &send, std::pair<T, size_t> &recv) {
        std::pair<T, size_t> tmp{send, irank()};
        return all_reduce(&tmp, &recv, MpiMaxLoc, 1);
    }

    /*
     * MINLOC CONVENIENCE FUNCTIONS
     */

    template<typename T>
    static bool all_minloc(const T &send, std::pair<T, size_t> &recv) {
        std::pair<T, size_t> tmp{send, irank()};
        return all_reduce(&tmp, &recv, MpiMinLoc, 1);
    }


    template<typename T>
    static bool all_to_all(const T *send, size_t nsend, T *recv, size_t nrecv) {
#ifdef HAVE_MPI
        return MPI_Alltoall((void *) send, snrw(nsend), mpi_type<T>(),
                            (void *) recv, snrw(nrecv), mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        ASSERT(nsend == nrecv);
        return all_reduce(send, recv, MpiMax, nsend);
#endif
    }

    template<typename T>
    static bool all_to_all(const std::vector<T> &send, std::vector<T> &recv) {
        return all_to_all(send.data(), 1ul, recv.data(), 1ul);
    }

private:
    template<typename T>
    static bool all_to_allv(
            const T *send, const defs::mpi_count *sendcounts, const defs::mpi_count *senddispls,
            T *recv, const defs::mpi_count *recvcounts, const defs::mpi_count *recvdispls) {
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
            const T *send, const defs::mpi_count& sendcount,
            T *recv, const defs::mpi_count *recvcounts, const defs::mpi_count *displs) {
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
        auto tmp_sendcounts = snrw(sendcounts);
        auto tmp_senddispls = snrw(senddispls);
        auto tmp_recvcounts = snrw(recvcounts);
        auto tmp_recvdispls = snrw(recvdispls);
        return all_to_allv(send, tmp_sendcounts.data(), tmp_senddispls.data(), recv,
                           tmp_recvcounts.data(), tmp_recvdispls.data());
    }

    template<typename T>
    static bool all_gatherv(
            const T *send, size_t sendcount,
            T *recv, const defs::inds &recvcounts, const defs::inds &recvdispls) {
        auto tmp_recvcounts = snrw(recvcounts);
        auto tmp_recvdispls = snrw(recvdispls);
        return all_gatherv(send, snrw(sendcount), recv, tmp_recvcounts.data(), tmp_recvdispls.data());
    }

    template<typename T>
    static bool all_gather(
            const T *send, size_t sendcount, T *recv, size_t recvcount) {
#ifdef HAVE_MPI
        return MPI_Allgather(
                (void *) send, snrw(sendcount), mpi_type<T>(), recv, snrw(recvcount), mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        return all_to_all(send, sendcount, recv, recvcount);
#endif
    }

    template<typename T>
    static bool all_gather(const T &send, std::vector<T> &recv) {
        if (recv.size() < nrank()) recv.resize(nrank());
        return all_gather(&send, 1ul, recv.data(), nrank());
    }

    template<typename T>
    static bool all_gatherv(const std::vector<T> &send, size_t sendcount, std::vector<T> &recv,
                            defs::inds& recvcounts, defs::inds& recvdispls) {
        recvcounts.resize(nrank());
        recvdispls.resize(nrank());
        all_gather(sendcount, recvcounts);
        counts_to_displs_consec(recvcounts, recvdispls);
        const size_t nrecv = recvdispls.back() + recvcounts.back();
        if (recv.size() < nrecv) recv.resize(nrecv);
        return all_gatherv(send.data(), sendcount, recv.data(), recvcounts, recvdispls);
    }

    template<typename T>
    static bool all_gatherv(const std::vector<T> &send, std::vector<T> &recv,
                            defs::inds& recvcounts, defs::inds& recvdispls) {
        return all_gatherv(send, send.size(), recv, recvcounts, recvdispls);
    }

    template<typename T>
    static bool all_gatherv(const std::vector<T> &send, size_t sendcount, std::vector<T> &recv) {
        defs::inds recvcounts, recvdispls;
        return all_gatherv(send, sendcount, recv, recvcounts, recvdispls);
    }

    template<typename T>
    static bool all_gatherv(const std::vector<T> &send, std::vector<T> &recv) {
        return all_gatherv(send, send.size(), recv);
    }

    static bool i_am(const size_t &i);

    static bool on_node_i_am(const size_t &i);

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
