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
#include "M7_lib/util/Convert.h"

/**
 * With MPI, we potentially have a typing issue. By default MPI libraries are typically
 * compiled with the displacement and count types set to 32-bit signed ints. This is
 * problematic whenever we want to communicate blocks of more than 2147483647 (2^32-1)
 * words (over 17GB worth of doubles). M7 uses unsigned longs (uint_t) throughout for
 * integer offsets, indices, and sizes - a mismatch this wrapper must gracefully resolve.
 *
 * The convert::safe_narrow methods will raise a runtime exception if a narrowing conversion
 * would result in a loss of information if the SAFE_NARROWING macro is defined. In normal
 * usage, such an overflow is fairly unlikely. But we do not leave it to chance that such a
 * large number of words will never be communicated in practice.
 *
 * count_t and countv_t are the scalar and vector typedefs for the MPI integer
 *
 * The cleanest way to achieve this aim is to expose only the overloads that use uint_t and
 * v_t<uint_t> (uintv_t) to define counts and displs. if sizeof(count_t)
 * is less than sizeof(uint_t), a narrowing conversion is required, but since this is at the
 * level of communication, and not inside a main loop, this minor overhead will more than
 * pay for itself in code clarity. Hence, **count_t and countv_t should never
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
static constexpr uint_t mpi_type_ind() { return ~0ul; }

template<>
constexpr uint_t mpi_type_ind<char>() { return 0; }

template<>
constexpr uint_t mpi_type_ind<short int>() { return 1; }

template<>
constexpr uint_t mpi_type_ind<int>() { return 2; }

template<>
constexpr uint_t mpi_type_ind<long int>() { return 3; }

template<>
constexpr uint_t mpi_type_ind<long long int>() { return 4; }

template<>
constexpr uint_t mpi_type_ind<unsigned char>() { return 5; }

template<>
constexpr uint_t mpi_type_ind<unsigned short int>() { return 6; }

template<>
constexpr uint_t mpi_type_ind<unsigned int>() { return 7; }

template<>
constexpr uint_t mpi_type_ind<unsigned long int>() { return 8; }

template<>
constexpr uint_t mpi_type_ind<unsigned long long int>() { return 9; }

template<>
constexpr uint_t mpi_type_ind<float>() { return 10; }

template<>
constexpr uint_t mpi_type_ind<double>() { return 11; }

template<>
constexpr uint_t mpi_type_ind<long double>() { return 12; }

template<>
constexpr uint_t mpi_type_ind<std::complex<float>>() { return 13; }

template<>
constexpr uint_t mpi_type_ind<std::complex<double>>() { return 14; }

template<>
constexpr uint_t mpi_type_ind<std::complex<long double>>() { return 15; }

template<>
constexpr uint_t mpi_type_ind<bool>() { return 16; }


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

/**
 * rank index in the world communicator
 */
extern uint_t g_irank;
/**
 * number of ranks in the world communicator
 */
extern uint_t g_nrank;
/**
 * name of this rank
 */
extern str_t g_processor_name;
/**
 * communicator among ranks on the same node, where a "node" is a group of ranks with access to the same main memory
 */
extern MPI_Comm g_node_comm;
/**
 * rank index within this node
 */
extern uint_t g_irank_on_node;
/**
 * number of ranks on this node
 */
extern uint_t g_nrank_on_node;
/**
 * list of world communicator rank indices of the node roots
 */
extern v_t<char> g_node_roots;
/**
 * world communicator rank index of the node root of this rank
 */
extern uint_t g_my_node_root_irank;
/**
 * todo: delete - point to point comms no longer used
 */
extern int g_p2p_tag;

/**
 * for performing multiple reduction operations simultaneously in a single communication (see reduction namespace)
 */
extern std::array<v_t<buf_t>, mpi_types.size()> g_send_reduction_buffers;
extern std::array<v_t<buf_t>, mpi_types.size()> g_recv_reduction_buffers;


namespace mpi {

    typedef int count_t;
    typedef v_t<count_t> countv_t;

    void setup_mpi_globals();

    static uint_t nrank() {
        return g_nrank;
    }

    static uint_t irank() {
        return g_irank;
    }

    static uint_t nrank_on_node() {
        return g_nrank_on_node;
    }

    static uint_t irank_on_node() {
        return g_irank_on_node;
    }

    static uint_t my_node_root_irank() {
        return g_my_node_root_irank;
    }

    static const str_t &processor_name() {
        return g_processor_name;
    }

    static MPI_Comm *node_communicator() {
        return &g_node_comm;
    }

    /*
     * NODE INDEXING
     */

    bool i_am(uint_t irank);

    bool on_node_i_am(uint_t irank);

    bool i_am_root();

    bool on_node_i_am_root();

    bool is_node_root(uint_t irank);

    void barrier();

    void barrier_on_node();

    static count_t snrw(uint_t i) {
        return convert::safe_narrow<count_t>(i);
    }

    static countv_t snrw(const uintv_t &v) {
        return convert::safe_narrow<count_t>(v);
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
    count_t evenly_shared_count(uint_t nitem_global, uint_t irank);
    count_t evenly_shared_count(uint_t nitem_global);
    countv_t evenly_shared_counts(uint_t nitem_global);

    uint_t evenly_shared_displ(uint_t nitem_global, uint_t irank);
    uint_t evenly_shared_displ(uint_t nitem_global);
    countv_t evenly_shared_displs(uint_t nitem_global);

    template<typename T>
    static void counts_to_displs_consec(const v_t<T> &counts, v_t<T> &displs) {
        ASSERT(counts.size() == displs.size())
        displs[0] = 0;
        for (uint_t i = 1ul; i < counts.size(); ++i) displs[i] = displs[i - 1] + counts[i - 1];
    }

    template<typename T>
    static v_t<T> counts_to_displs_consec(const v_t<T> &counts) {
        auto displs = counts;
        counts_to_displs_consec(counts, displs);
        return displs;
    }

    template<typename T>
    static bool reduce(const T *send, T *recv, MpiOp op, uint_t ndata = 1, uint_t iroot = 0) {
        return MPI_Reduce(send, recv, snrw(ndata), mpi_type<T>(), op_map[op], snrw(iroot), MPI_COMM_WORLD) ==
               MPI_SUCCESS;
    }

    template<typename T>
    static bool all_reduce(const T *send, T *recv, MpiOp op, uint_t ndata = 1) {
        static_assert(mpi_type_ind<T>() != ~0ul, "Not a valid MPI type");
        return MPI_Allreduce(send, recv, snrw(ndata), mpi_type<T>(), op_map[op], MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool all_reduce(const std::pair<T, uint_t> *send,
                           std::pair<T, uint_t> *recv,
                           MpiPairOp op, uint_t ndata = 1) {
        static_assert(mpi_type_ind<T>() != ~0ul, "Not a valid MPI type");
        v_t<std::pair<T, count_t>> tmp_send;
        tmp_send.reserve(ndata);
        for (uint_t idata = 0ul; idata < ndata; ++idata)
            tmp_send.push_back({send[idata].first, snrw(send[idata].second)});
        v_t<std::pair<T, count_t>> tmp_recv;
        tmp_recv.reserve(ndata);
        auto res = MPI_Allreduce(tmp_send.data(), tmp_recv.data(), ndata, mpi_pair_type<T>(),
                                 pair_op_map[op], MPI_COMM_WORLD) == MPI_SUCCESS;
        for (uint_t idata = 0ul; idata < ndata; ++idata)
            recv[idata] = std::pair<T, uint_t>{tmp_recv[idata].first, tmp_recv[idata].second};
        return res;
    }


    template<typename T>
    static bool all_reduce(const v_t<std::pair<T, uint_t>> &send,
                           v_t<std::pair<T, uint_t>> &recv,
                           MpiPairOp op) {
        ASSERT(send.size() == recv.size());
        return all_reduce(send.data(), recv.data(), op, send.size());
    }

    /*
     * MAX REDUCE CONVENIENCE METHODS
     */
    template<typename T>
    static bool max(const T *send, T *recv, uint_t ndata = 1, uint_t iroot = 0) {
        return reduce(send, recv, MpiMax, ndata, iroot);
    }

    template<typename T>
    static bool all_max(const T *send, T *recv, uint_t ndata = 1) {
        return all_reduce(send, recv, MpiMax, ndata);
    }

    template<typename T>
    static T max(const T *send, uint_t iroot = 0) {
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
    static bool min(const T *send, T *recv, uint_t ndata = 1, uint_t iroot = 0) {
        return reduce(send, recv, MpiMin, ndata, iroot);
    }

    template<typename T>
    static bool all_min(const T *send, T *recv, uint_t ndata = 1) {
        return all_reduce(send, recv, MpiMin, ndata);
    }

    template<typename T>
    static T min(const T *send, uint_t iroot = 0) {
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
    static bool sum(const T *send, T *recv, uint_t ndata = 1, uint_t iroot = 0) {
        return reduce(send, recv, MpiSum, ndata, iroot);
    }

    template<typename T>
    static bool all_sum(const T *send, T *recv, uint_t ndata = 1) {
        return all_reduce(send, recv, MpiSum, ndata);
    }

    template<typename T>
    static T sum(const T *send, uint_t iroot = 0) {
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
    static bool land(const T *send, T *recv, uint_t ndata = 1, uint_t iroot = 0) {
        return reduce(send, recv, MpiLand, ndata, iroot);
    }

    template<typename T>
    static bool all_land(const T *send, T *recv, uint_t ndata = 1) {
        return all_reduce(send, recv, MpiLand, ndata);
    }

    template<typename T>
    static T land(const T *send, uint_t iroot = 0) {
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
    static bool lor(const T *send, T *recv, uint_t ndata = 1, uint_t iroot = 0) {
        return reduce(send, recv, MpiLor, ndata, iroot);
    }

    template<typename T>
    static bool all_lor(const T *send, T *recv, uint_t ndata = 1) {
        return all_reduce(send, recv, MpiLor, ndata);
    }

    template<typename T>
    static T lor(const T *send, uint_t iroot = 0) {
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
    static bool send(const T *data, uint_t ndata, uint_t irank_dst, int tag) {
        return MPI_Send(reinterpret_cast<const void *>(data), snrw(ndata), mpi_type<T>(), snrw(irank_dst), tag,
                        MPI_COMM_WORLD);
    }

    template<typename T>
    static bool recv(T *data, uint_t ndata, uint_t irank_src, int tag) {
        return MPI_Recv(reinterpret_cast<void *>(data), snrw(ndata), mpi_type<T>(), snrw(irank_src), tag,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE);
    }


    /*
     * BCAST
     */
    template<typename T>
    static bool bcast(T *data, uint_t ndata = 1, uint_t iroot = 0) {
        return MPI_Bcast(reinterpret_cast<void *>(data), snrw(ndata), mpi_type<T>(), snrw(iroot), MPI_COMM_WORLD) ==
               MPI_SUCCESS;
    }

    template<typename T>
    static bool bcast(T &data, uint_t iroot = 0) {
        return bcast(&data, 1, iroot);
    }

    template<typename T>
    static bool bcast(v_t<T> &data, uint_t iroot = 0) {
        uint_t size = data.size();
        mpi::bcast(size, iroot);
        if (!i_am(iroot)) data.resize(size);
        auto data_ptr = reinterpret_cast<void *>(data.data());
        return MPI_Bcast(data_ptr, snrw(size), mpi_type<T>(), snrw(iroot), MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    static bool bcast(str_t &data, uint_t iroot = 0) {
        uint_t nchar = data.size();
        mpi::bcast(nchar, iroot);
        data.resize(nchar);
        return mpi::bcast(const_cast<char *>(data.c_str()), nchar);
    }

    /*
     * MAXLOC CONVENIENCE FUNCTIONS
     */
    template<typename T>
    static bool all_maxloc(const T *send, std::pair<T, uint_t> *recv, uint_t ndata) {
        v_t<std::pair<T, uint_t>> tmp;
        tmp.reserve(ndata);
        for (uint_t idata = 0ul; idata < ndata; ++idata) tmp.push_back({send[idata], irank()});
        return all_reduce(tmp.data(), recv, MpiMaxLoc, ndata);
    }

    template<typename T>
    static bool all_maxloc(const T &send, std::pair<T, uint_t> &recv) {
        return all_maxloc(&send, &recv, 1);
    }

    template<typename T>
    static bool all_maxloc(const v_t<T> &send, v_t<std::pair<T, uint_t>> &recv) {
        ASSERT(send.size() == recv.size());
        return all_maxloc(send.data(), recv.data(), send.size());
    }

    /*
     * MINLOC CONVENIENCE FUNCTIONS
     */
    template<typename T>
    static bool all_minloc(const T *send, std::pair<T, uint_t> *recv, uint_t ndata) {
        v_t<std::pair<T, uint_t>> tmp;
        tmp.reserve(ndata);
        for (uint_t idata = 0ul; idata < ndata; ++idata) tmp.push_back({send[idata], irank()});
        return all_reduce(tmp.data(), recv, MpiMinLoc, ndata);
    }

    template<typename T>
    static bool all_minloc(const T &send, std::pair<T, uint_t> &recv) {
        return all_minloc(&send, &recv, 1);
    }

    template<typename T>
    static bool all_minloc(const v_t<T> &send, v_t<std::pair<T, uint_t>> &recv) {
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
            const T *send, uint_t sendcount, T *recv, uint_t recvcount, uint_t iroot) {
        auto send_ptr = reinterpret_cast<void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Gather(send_ptr, snrw(sendcount), mpi_type<T>(), recv_ptr,
                          snrw(recvcount), mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool all_gather(
            const T *send, uint_t sendcount, T *recv, uint_t recvcount) {
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
    static bool all_gather(const T &send, v_t<T> &recv) {
        recv.resize(nrank());
        return all_gather(send, recv.data());
    }

    template<typename T>
    static bool gatherv(
            const T *send, const count_t &sendcount,
            T *recv, const count_t *recvcounts, const count_t *recvdispls, uint_t iroot) {
        auto send_ptr = reinterpret_cast<const void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Gatherv(
                send_ptr, sendcount, mpi_type<T>(),
                recv_ptr, recvcounts, recvdispls, mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool all_gatherv(
            const T *send, const count_t &sendcount,
            T *recv, const count_t *recvcounts, const count_t *recvdispls) {
        auto send_ptr = reinterpret_cast<const void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Allgatherv(
                send_ptr, sendcount, mpi_type<T>(),
                recv_ptr, recvcounts, recvdispls, mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool gatherv(const T *send, uint_t sendcount, T *recv,
                        const uintv_t &recvcounts, const uintv_t &recvdispls, uint_t iroot) {
        auto tmp_recvcounts = snrw(recvcounts);
        auto tmp_recvdispls = snrw(recvdispls);
        return gatherv(send, snrw(sendcount), recv, tmp_recvcounts.data(), tmp_recvdispls.data(), iroot);
    }

    template<typename T>
    static bool all_gatherv(
            const T *send, uint_t sendcount,
            T *recv, const uintv_t &recvcounts, const uintv_t &recvdispls) {
        auto tmp_recvcounts = snrw(recvcounts);
        auto tmp_recvdispls = snrw(recvdispls);
        return all_gatherv(send, snrw(sendcount), recv, tmp_recvcounts.data(), tmp_recvdispls.data());
    }

    template<typename T>
    static bool all_gatherv(const v_t<T>& send, v_t<T>& recv) {
        uintv_t counts;
        all_gather(send.size(), counts);
        const auto displs = counts_to_displs_consec(counts);
        const auto n = displs.back() + counts.back();
        recv.resize(n);
        return all_gatherv(send.data(), send.size(), recv.data(), counts, displs);
    }

    template<typename T>
    static v_t<T> all_gathered(const T& send) {
        v_t<T> out(mpi::nrank());
        all_gather(send, out);
        return out;
    }

    template<typename T>
    static v_t<T> all_gatheredv(const v_t<T>& send) {
        v_t<T> out;
        all_gatherv(send, out);
        return out;
    }

    /**
     * @param cond
     *  true if this rank should be included in iranks at output
     * @param iranks
     *  ordered vector of rank indices that pass the filter condition
     */
    void filter(bool cond, uintv_t& iranks);

    uintv_t filter(bool cond);

    /*
     * SCATTER
     */

    template<typename T>
    static bool scatter(
            const T *send, uint_t sendcount, T *recv, uint_t recvcount, uint_t iroot) {
        auto send_ptr = reinterpret_cast<void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Scatter(send_ptr, snrw(sendcount), mpi_type<T>(), recv_ptr,
                           snrw(recvcount), mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool scatterv(
            const T *send, const count_t *sendcounts, const count_t *senddispls,
            T *recv, const count_t &recvcount, uint_t iroot) {
        auto send_ptr = reinterpret_cast<const void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Scatterv(
                send_ptr, sendcounts, senddispls, mpi_type<T>(),
                recv_ptr, recvcount, mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool scatterv(const T *send, const uintv_t sendcounts, const uintv_t &senddispls,
                         T *recv, uint_t recvcount, uint_t iroot) {
        auto tmp_sendcounts = snrw(sendcounts);
        auto tmp_senddispls = snrw(senddispls);
        return gatherv(send, tmp_sendcounts, tmp_senddispls, recv, snrw(recvcount), iroot);
    }


    /*
     * ALL-TO-ALL
     */

    template<typename T>
    static bool all_to_all(const T *send, uint_t nsend, T *recv, uint_t nrecv) {
        auto send_ptr = reinterpret_cast<const void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Alltoall(send_ptr, snrw(nsend), mpi_type<T>(),
                            recv_ptr, snrw(nrecv), mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool all_to_all(const v_t<T> &send, v_t<T> &recv) {
        return all_to_all(send.data(), 1ul, recv.data(), 1ul);
    }

    template<typename T>
    static bool all_to_allv(
            const T *send, const count_t *sendcounts, const count_t *senddispls,
            T *recv, const count_t *recvcounts, const count_t *recvdispls) {
        auto send_ptr = reinterpret_cast<const void *>(send);
        auto recv_ptr = reinterpret_cast<void *>(recv);
        return MPI_Alltoallv(
                send_ptr, sendcounts, senddispls, mpi_type<T>(),
                recv_ptr, recvcounts, recvdispls, mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
    }

    template<typename T>
    static bool all_to_allv(
            const T *send, const uintv_t &sendcounts, const uintv_t &senddispls,
            T *recv, const uintv_t &recvcounts, const uintv_t &recvdispls) {
        auto tmp_sendcounts = snrw(sendcounts);
        auto tmp_senddispls = snrw(senddispls);
        auto tmp_recvcounts = snrw(recvcounts);
        auto tmp_recvdispls = snrw(recvdispls);
        return all_to_allv(send, tmp_sendcounts.data(), tmp_senddispls.data(), recv,
                           tmp_recvcounts.data(), tmp_recvdispls.data());
    }

    bool initialized();

    bool finalized();

    void initialize(int *argc = nullptr, char ***argv = nullptr);

    void finalize();

    /**
     * call on any rank to terminate the entire communicator
     * @param message
     *  string to output to error log
     */
    void abort_(str_t message);

    /**
     * call on all ranks to terminate the entire communicator
     * @param message
     *  string to output to error log
     */
    void abort(str_t message);

    /**
     * debugging: each rank waits till the one before has finished before starting
     * @param str
     *  string to output to stdout
     */
    void blocking_print(const str_t &str);

}

#endif //M7_MPIWRAPPER_H
