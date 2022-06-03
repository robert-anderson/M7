//
// Created by Robert John Anderson on 2020-02-07.
//

#ifndef M7_MPIWRAPPER_H
#define M7_MPIWRAPPER_H

#include <array>

#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif

#include <mpi.h>
// we don't use the C++ MPI bindings

#include <iostream>
#include <cstring>

#include <M7_lib/defs.h>
#include <M7_lib/util/utils.h>

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
 * The cleanest way to achieve this aim is to expose only the overloads that use size_t and
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

static const std::array<MPI_Datatype, 17> mpi_types =
        {MPI_CHAR, MPI_SHORT, MPI_INT, MPI_LONG, MPI_LONG_LONG_INT,
         MPI_UNSIGNED_CHAR, MPI_UNSIGNED_SHORT, MPI_UNSIGNED, MPI_UNSIGNED_LONG,
         MPI_UNSIGNED_LONG_LONG, MPI_FLOAT, MPI_DOUBLE, MPI_LONG_DOUBLE,
         MPI_COMPLEX, MPI_DOUBLE_COMPLEX, MPI_CXX_LONG_DOUBLE_COMPLEX, MPI_CXX_BOOL};

template<typename T>
static constexpr size_t mpi_type_ind() { return ~0ul; }

template<>
constexpr size_t mpi_type_ind<char>() { return 0; }

template<>
constexpr size_t mpi_type_ind<short int>() { return 1; }

template<>
constexpr size_t mpi_type_ind<int>() { return 2; }

template<>
constexpr size_t mpi_type_ind<long int>() { return 3; }

template<>
constexpr size_t mpi_type_ind<long long int>() { return 4; }

template<>
constexpr size_t mpi_type_ind<unsigned char>() { return 5; }

template<>
constexpr size_t mpi_type_ind<unsigned short int>() { return 6; }

template<>
constexpr size_t mpi_type_ind<unsigned int>() { return 7; }

template<>
constexpr size_t mpi_type_ind<unsigned long int>() { return 8; }

template<>
constexpr size_t mpi_type_ind<unsigned long long int>() { return 9; }

template<>
constexpr size_t mpi_type_ind<float>() { return 10; }

template<>
constexpr size_t mpi_type_ind<double>() { return 11; }

template<>
constexpr size_t mpi_type_ind<long double>() { return 12; }

template<>
constexpr size_t mpi_type_ind<std::complex<float>>() { return 13; }

template<>
constexpr size_t mpi_type_ind<std::complex<double>>() { return 14; }

template<>
constexpr size_t mpi_type_ind<std::complex<long double>>() { return 15; }

template<>
constexpr size_t mpi_type_ind<bool>() { return 16; }


template<typename T>
constexpr bool mpi_supported() { return mpi_type_ind<T>() != ~0ul; }


template<typename T>
static MPI_Datatype mpi_type() { return mpi_types[mpi_type_ind<T>()]; }

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

extern size_t g_irank;
extern size_t g_nrank;
extern std::string g_processor_name;
extern MPI_Comm g_node_comm;
extern size_t g_irank_on_node;
extern size_t g_nrank_on_node;
extern int g_p2p_tag;

namespace mpi {

    void setup_mpi_globals();

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
        return &g_node_comm;
    }

    void barrier();

    void barrier_on_node();

    static defs::mpi_count snrw(const size_t &i) {
        return utils::safe_narrow<defs::mpi_count>(i);
    }

    static defs::mpi_counts snrw(const std::vector<size_t> &v) {
        return utils::safe_narrow<defs::mpi_count>(v);
    }

    /**
     *
     * e.g. nrank = 6, nitem = 16
     *  remainder = 16 % 6 = 4
     *  count = nitem/nrank + (irank < remainder) =
     *         3  3  3  3  2  2
     *
     *  number of tall stacks before irank: ntb = min(irank, remainder)
     *  number of short stacks before irank: nsb = irank - ntb
     *  displ = ntb * (ss + 1) + nsb * ss
     *        = ntb + irank * ss =
     *         0  3  6  9 12 14
     *
     *  rank   0  1  2  3  4  5
     *         x  x  x  x  x  x
     *         x  x  x  x  x  x
     *         x  x  x  x
     * @param nitem_global
     * @return
     */
    defs::mpi_count evenly_shared_count(size_t nitem_global, size_t irank);
    defs::mpi_count evenly_shared_count(size_t nitem_global);
    defs::mpi_counts evenly_shared_counts(size_t nitem_global);

    size_t evenly_shared_displ(size_t nitem_global, size_t irank);
    size_t evenly_shared_displ(size_t nitem_global);
    defs::mpi_counts evenly_shared_displs(size_t nitem_global);

    template<typename T>
    static void counts_to_displs_consec(const std::vector<T> &counts, std::vector<T> &displs) {
        ASSERT(counts.size() == displs.size())
        displs[0] = 0;
        for (size_t i = 1ul; i < counts.size(); ++i) displs[i] = displs[i - 1] + counts[i - 1];
    }

    template<typename T>
    static std::vector<T> counts_to_displs_consec(const std::vector<T> &counts) {
        auto displs = counts;
        counts_to_displs_consec(counts, displs);
        return displs;
    }

    template<typename T>
    static bool reduce(const T *send, T *recv, MpiOp op, size_t ndata = 1, size_t iroot = 0) {
        return MPI_Reduce(send, recv, snrw(ndata), mpi_type<T>(), op_map[op], snrw(iroot), MPI_COMM_WORLD) ==
               MPI_SUCCESS;
    }

    template<typename T>
    static bool all_reduce(const T *send, T *recv, MpiOp op, size_t ndata = 1) {
        static_assert(mpi_type_ind<T>() != ~0ul, "Not a valid MPI type");
        return MPI_Allreduce(send, recv, snrw(ndata), mpi_type<T>(), op_map[op], MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool all_reduce(const std::pair<T, size_t> *send,
                           std::pair<T, size_t> *recv,
                           MpiPairOp op, size_t ndata = 1) {
        static_assert(mpi_type_ind<T>() != ~0ul, "Not a valid MPI type");
        std::vector<std::pair<T, defs::mpi_count>> tmp_send;
        tmp_send.reserve(ndata);
        for (size_t idata = 0ul; idata < ndata; ++idata)
            tmp_send.push_back({send[idata].first, snrw(send[idata].second)});
        std::vector<std::pair<T, defs::mpi_count>> tmp_recv;
        tmp_recv.reserve(ndata);
        auto res = MPI_Allreduce(tmp_send.data(), tmp_recv.data(), ndata, mpi_pair_type<T>(),
                                 pair_op_map[op], MPI_COMM_WORLD) == MPI_SUCCESS;
        for (size_t idata = 0ul; idata < ndata; ++idata)
            recv[idata] = std::pair<T, size_t>{tmp_recv[idata].first, tmp_recv[idata].second};
        return res;
    }


    template<typename T>
    static bool all_reduce(const std::vector<std::pair<T, size_t>> &send,
                           std::vector<std::pair<T, size_t>> &recv,
                           MpiPairOp op) {
        ASSERT(send.size() == recv.size());
        return all_reduce(send.data(), recv.data(), op, send.size());
    }

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
    static T all_land(T send) {
        T res;
        T cpy = send;
        all_land(&cpy, &res);
        return res;
    }

    /*
     * LOGICAL-OR REDUCE CONVENIENCE FUNCTIONS
     */

    template<typename T>
    static bool lor(const T *send, T *recv, size_t ndata = 1, size_t iroot = 0) {
        return reduce(send, recv, MpiLor, ndata, iroot);
    }

    template<typename T>
    static bool all_lor(const T *send, T *recv, size_t ndata = 1) {
        return all_reduce(send, recv, MpiLor, ndata);
    }

    template<typename T>
    static T lor(const T *send, size_t iroot = 0) {
        T recv;
        lor(send, &recv, 1, iroot);
        return recv;
    }

    template<typename T>
    static T all_lor(const T *send) {
        T recv;
        all_lor(send, &recv, 1);
        return recv;
    }

    template<typename T>
    static T all_lor(T send) {
        T res;
        all_lor(&send, &res);
        return res;
    }


    /*
     * P2P send / recv
     * We require a globally-unique tag for each P2P send/recv pair in the whole program.
     */
    static int new_p2p_tag() {
        return g_p2p_tag++;
    }

    template<typename T>
    static bool send(const T *data, size_t ndata, size_t irank_dst, int tag) {
        return MPI_Send(reinterpret_cast<const void *>(data), snrw(ndata), mpi_type<T>(), snrw(irank_dst), tag,
                        MPI_COMM_WORLD);
    }

    template<typename T>
    static bool recv(T *data, size_t ndata, size_t irank_src, int tag) {
        return MPI_Recv(reinterpret_cast<void *>(data), snrw(ndata), mpi_type<T>(), snrw(irank_src), tag,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE);
    }


    /*
     * BCAST
     */
    template<typename T>
    static bool bcast(T *data, size_t ndata = 1, size_t iroot = 0) {
        return MPI_Bcast(reinterpret_cast<void *>(data), snrw(ndata), mpi_type<T>(), snrw(iroot), MPI_COMM_WORLD) ==
               MPI_SUCCESS;
    }

    template<typename T>
    static bool bcast(T &data, size_t iroot = 0) {
        return bcast(&data, 1, iroot);
    }

    template<typename T>
    static bool bcast(std::vector<T> &data, size_t ndata = 0, size_t iroot = 0) {
        auto data_ptr = reinterpret_cast<void *>(data.data());
        if (!ndata) ndata = data.size();
        data.resize(ndata);
        return MPI_Bcast(data_ptr, snrw(ndata), mpi_type<T>(), snrw(iroot), MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    static bool bcast(std::string &data, size_t iroot = 0) {
        size_t nchar = data.size();
        mpi::bcast(nchar, iroot);
        data.resize(nchar);
        return mpi::bcast(const_cast<char *>(data.c_str()), nchar);
    }

    /*
     * MAXLOC CONVENIENCE FUNCTIONS
     */
    template<typename T>
    static bool all_maxloc(const T *send, std::pair<T, size_t> *recv, size_t ndata) {
        std::vector<std::pair<T, size_t>> tmp;
        tmp.reserve(ndata);
        for (size_t idata = 0ul; idata < ndata; ++idata) tmp.push_back({send[idata], irank()});
        return all_reduce(tmp.data(), recv, MpiMaxLoc, ndata);
    }

    template<typename T>
    static bool all_maxloc(const T &send, std::pair<T, size_t> &recv) {
        return all_maxloc(&send, &recv, 1);
    }

    template<typename T>
    static bool all_maxloc(const std::vector<T> &send, std::vector<std::pair<T, size_t>> &recv) {
        ASSERT(send.size() == recv.size());
        return all_maxloc(send.data(), recv.data(), send.size());
    }

    /*
     * MINLOC CONVENIENCE FUNCTIONS
     */
    template<typename T>
    static bool all_minloc(const T *send, std::pair<T, size_t> *recv, size_t ndata) {
        std::vector<std::pair<T, size_t>> tmp;
        tmp.reserve(ndata);
        for (size_t idata = 0ul; idata < ndata; ++idata) tmp.push_back({send[idata], irank()});
        return all_reduce(tmp.data(), recv, MpiMinLoc, ndata);
    }

    template<typename T>
    static bool all_minloc(const T &send, std::pair<T, size_t> &recv) {
        return all_minloc(&send, &recv, 1);
    }

    template<typename T>
    static bool all_minloc(const std::vector<T> &send, std::vector<std::pair<T, size_t>> &recv) {
        ASSERT(send.size() == recv.size());
        return all_minloc(send.data(), recv.data(), send.size());
    }

    /*
     * COLLECTIVE CALLS
     */

    /*
     * GATHER
     */

    template<typename T>
    static bool gather(
            const T *send, size_t sendcount, T *recv, size_t recvcount, const size_t &iroot) {
        auto send_ptr = reinterpret_cast<void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Gather(send_ptr, snrw(sendcount), mpi_type<T>(), recv_ptr,
                          snrw(recvcount), mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool all_gather(
            const T *send, size_t sendcount, T *recv, size_t recvcount) {
        auto send_ptr = reinterpret_cast<const void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Allgather(
                send_ptr, snrw(sendcount), mpi_type<T>(), recv_ptr, snrw(recvcount), mpi_type<T>(), MPI_COMM_WORLD) ==
               MPI_SUCCESS;
    }

    template<typename T>
    static bool all_gather(const T &send, T *recv) {
        return all_gather(&send, 1, recv, 1);
    }

    template<typename T>
    static bool all_gather(const T &send, std::vector<T> &recv) {
        recv.resize(nrank());
        return all_gather(send, recv.data());
    }

    template<typename T>
    static bool gatherv(
            const T *send, const defs::mpi_count &sendcount,
            T *recv, const defs::mpi_count *recvcounts, const defs::mpi_count *recvdispls, const size_t &iroot) {
        auto send_ptr = reinterpret_cast<const void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Gatherv(
                send_ptr, sendcount, mpi_type<T>(),
                recv_ptr, recvcounts, recvdispls, mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool all_gatherv(
            const T *send, const defs::mpi_count &sendcount,
            T *recv, const defs::mpi_count *recvcounts, const defs::mpi_count *recvdispls) {
        auto send_ptr = reinterpret_cast<const void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Allgatherv(
                send_ptr, sendcount, mpi_type<T>(),
                recv_ptr, recvcounts, recvdispls, mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool gatherv(const T *send, size_t sendcount, T *recv,
                        const defs::inds &recvcounts, const defs::inds &recvdispls, const size_t &iroot) {
        auto tmp_recvcounts = snrw(recvcounts);
        auto tmp_recvdispls = snrw(recvdispls);
        return gatherv(send, snrw(sendcount), recv, tmp_recvcounts.data(), tmp_recvdispls.data(), iroot);
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
    static std::vector<T> all_gathered(const T &send) {
        std::vector<T> out(mpi::nrank());
        all_gather(send, out);
        return out;
    }

    /*
     * SCATTER
     */

    template<typename T>
    static bool scatter(
            const T *send, size_t sendcount, T *recv, size_t recvcount, const size_t &iroot) {
        auto send_ptr = reinterpret_cast<void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Scatter(send_ptr, snrw(sendcount), mpi_type<T>(), recv_ptr,
                           snrw(recvcount), mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool scatterv(
            const T *send, const defs::mpi_count *sendcounts, const defs::mpi_count *senddispls,
            T *recv, const defs::mpi_count &recvcount, const size_t &iroot) {
        auto send_ptr = reinterpret_cast<const void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Scatterv(
                send_ptr, sendcounts, senddispls, mpi_type<T>(),
                recv_ptr, recvcount, mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool scatterv(const T *send, const defs::inds sendcounts, const defs::inds &senddispls,
                         T *recv, const size_t &recvcount, const size_t &iroot) {
        auto tmp_sendcounts = snrw(sendcounts);
        auto tmp_senddispls = snrw(senddispls);
        return gatherv(send, tmp_sendcounts, tmp_senddispls, recv, snrw(recvcount), iroot);
    }


    /*
     * ALL-TO-ALL
     */

    template<typename T>
    static bool all_to_all(const T *send, size_t nsend, T *recv, size_t nrecv) {
        auto send_ptr = reinterpret_cast<const void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Alltoall(send_ptr, snrw(nsend), mpi_type<T>(),
                            recv_ptr, snrw(nrecv), mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool all_to_all(const std::vector<T> &send, std::vector<T> &recv) {
        return all_to_all(send.data(), 1ul, recv.data(), 1ul);
    }

    template<typename T>
    static bool all_to_allv(
            const T *send, const defs::mpi_count *sendcounts, const defs::mpi_count *senddispls,
            T *recv, const defs::mpi_count *recvcounts, const defs::mpi_count *recvdispls) {
        auto send_ptr = reinterpret_cast<const void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Alltoallv(
                send_ptr, sendcounts, senddispls, mpi_type<T>(),
                recv_ptr, recvcounts, recvdispls, mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
    }

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

    /*
     * NODE INDEXING
     */

    bool i_am(const size_t &i);

    bool on_node_i_am(const size_t &i);

    bool i_am_root();

    bool on_node_i_am_root();

    bool initialized();

    bool finalized();

    void initialize(int *argc = nullptr, char ***argv = nullptr);

    void finalize();

    /**
     * call on any rank to terminate the entire communicator
     * @param message
     *  string to output to error log
     */
    void abort_(std::string message);

    /**
     * call on all ranks to terminate the entire communicator
     * @param message
     *  string to output to error log
     */
    void abort(std::string message);

    /**
     * debugging: each rank waits till the one before has finished before starting
     * @param str
     *  string to output to stdout
     */
    void blocking_print(const std::string &str);

}

#endif //M7_MPIWRAPPER_H
