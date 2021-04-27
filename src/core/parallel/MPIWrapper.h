//
// Created by Robert John Anderson on 2020-02-07.
//

#ifndef M7_MPIWRAPPER_H
#define M7_MPIWRAPPER_H

#include <array>

#ifdef ENABLE_MPI
#define OMPI_SKIP_MPICXX
#include <mpi.h>
// we don't use the C++ bindings
#endif

#include <iostream>
#include <cstring>
#include <src/defs.h>
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
#ifdef ENABLE_MPI


static const std::array<MPI_Datatype, 17> mpi_types =
        {MPI_CHAR, MPI_SHORT, MPI_INT, MPI_LONG, MPI_LONG_LONG_INT,
         MPI_UNSIGNED_CHAR, MPI_UNSIGNED_SHORT, MPI_UNSIGNED, MPI_UNSIGNED_LONG,
         MPI_UNSIGNED_LONG_LONG, MPI_FLOAT, MPI_DOUBLE, MPI_LONG_DOUBLE,
         MPI_COMPLEX, MPI_DOUBLE_COMPLEX, MPI_CXX_LONG_DOUBLE_COMPLEX, MPI_CXX_BOOL};

template<typename T>
static constexpr size_t mpi_type_ind() { return ~0ul; }

template<> constexpr size_t mpi_type_ind<char>() { return 0;}
template<> constexpr size_t mpi_type_ind<short int>() { return 1;}
template<> constexpr size_t mpi_type_ind<int>() { return 2; }
template<> constexpr size_t mpi_type_ind<long int>() { return 3; }
template<> constexpr size_t mpi_type_ind<long long int>() { return 4; }
template<> constexpr size_t mpi_type_ind<unsigned char>() { return 5; }
template<> constexpr size_t mpi_type_ind<unsigned short int>() { return 6; }
template<> constexpr size_t mpi_type_ind<unsigned int>() { return 7; }
template<> constexpr size_t mpi_type_ind<unsigned long int>() { return 8; }
template<> constexpr size_t mpi_type_ind<unsigned long long int>() { return 9; }
template<> constexpr size_t mpi_type_ind<float>() { return 10; }
template<> constexpr size_t mpi_type_ind<double>() { return 11; }
template<> constexpr size_t mpi_type_ind<long double>() { return 12; }
template<> constexpr size_t mpi_type_ind<std::complex<float>>() { return 13; }
template<> constexpr size_t mpi_type_ind<std::complex<double>>() { return 14; }
template<> constexpr size_t mpi_type_ind<std::complex<long double>>() { return 15; }
template<> constexpr size_t mpi_type_ind<bool>() { return 16; }


template<typename T> constexpr bool mpi_supported() { return mpi_type_ind<T>()!=~0ul; }


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

template<> MPI_Datatype mpi_pair_type<float>() { return MPI_FLOAT_INT; }
template<> MPI_Datatype mpi_pair_type<long>() { return MPI_LONG_INT; }
template<> MPI_Datatype mpi_pair_type<double>() { return MPI_DOUBLE_INT; }
template<> MPI_Datatype mpi_pair_type<short>() { return MPI_SHORT_INT; }
template<> MPI_Datatype mpi_pair_type<int>() { return MPI_2INT; }
template<> MPI_Datatype mpi_pair_type<long double>() { return MPI_LONG_DOUBLE_INT; }

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
#ifdef ENABLE_MPI
extern MPI_Comm g_node_comm;
#endif
extern size_t g_irank_on_node;
extern size_t g_nrank_on_node;
extern int g_p2p_tag;

struct mpi {

    static void setup_mpi_globals();

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
#ifdef ENABLE_MPI
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
#ifdef ENABLE_MPI
        return MPI_Reduce(send, recv, snrw(ndata), mpi_type<T>(), op_map[op], snrw(iroot), MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        std::memcpy(recv, send, sizeof(T)*ndata);
        return true;
#endif
    }

    template<typename T>
    static bool all_reduce(const T *send, T *recv, MpiOp op, size_t ndata = 1) {
#ifdef ENABLE_MPI
        static_assert(mpi_type_ind<T>() != ~0ul, "Not a valid MPI type");
        return MPI_Allreduce(send, recv, snrw(ndata), mpi_type<T>(), op_map[op], MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        std::memcpy(recv, send, sizeof(T)*ndata);
        return true;
#endif
    }

    template<typename T>
    static bool all_reduce(const std::pair<T, size_t>* send,
                           std::pair<T, size_t>* recv,
                           MpiPairOp op, size_t ndata=1) {
#ifdef ENABLE_MPI
        static_assert(mpi_type_ind<T>() != ~0ul, "Not a valid MPI type");
        std::vector<std::pair<T, defs::mpi_count>> tmp_send;
        tmp_send.reserve(ndata);
        for (size_t idata=0ul; idata<ndata; ++idata)
            tmp_send.push_back({send[idata].first, snrw(send[idata].second)});
        std::vector<std::pair<T, defs::mpi_count>> tmp_recv;
        tmp_recv.reserve(ndata);
        auto res = MPI_Allreduce(tmp_send.data(), tmp_recv.data(), ndata, mpi_pair_type<T>(),
                pair_op_map[op], MPI_COMM_WORLD) == MPI_SUCCESS;
        for (size_t idata=0ul; idata<ndata; ++idata)
            recv[idata] = std::pair<T, size_t>{tmp_recv[idata].first, tmp_recv[idata].second};
        return res;
#else
        std::memcpy(recv.data(), send.data(), sizeof(std::pair<T, size_t>)*send.size());
        return true;
#endif
    }


    template<typename T>
    static bool all_reduce(const std::vector<std::pair<T, size_t>> &send,
                           std::vector<std::pair<T, size_t>> &recv,
                           MpiPairOp op) {
        ASSERT(send.size() == recv.size());
        return all_reduce(send.data(), recv.data(), op, send.size());
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
#ifdef ENABLE_MPI
        return MPI_Send((void*) data, snrw(ndata), mpi_type<T>(), snrw(irank_dst), tag, MPI_COMM_WORLD);
#endif
    }

    template<typename T>
    static bool recv(T *data, size_t ndata, size_t irank_src, int tag) {
#ifdef ENABLE_MPI
        return MPI_Recv((void*) data, snrw(ndata), mpi_type<T>(), snrw(irank_src), tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
    }


    /*
     * BCAST
     */
    template<typename T>
    static bool bcast(T *data, size_t ndata = 1, size_t iroot = 0) {
#ifdef ENABLE_MPI
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
#ifdef ENABLE_MPI
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
    static bool all_maxloc(const T *send, std::pair<T, size_t>* recv, size_t ndata) {
        std::vector<std::pair<T, size_t>> tmp;
        tmp.reserve(ndata);
        for (size_t idata=0ul; idata<ndata; ++idata) tmp.push_back({send[idata], irank()});
        return all_reduce(tmp.data(), recv, MpiMaxLoc, ndata);
    }

    template<typename T>
    static bool all_maxloc(const T &send, std::pair<T, size_t> &recv) {
        return all_maxloc(&send, &recv, 1);
    }

    template<typename T>
    static bool all_maxloc(const std::vector<T> &send, std::vector<std::pair<T, size_t>> &recv) {
        ASSERT(send.size()==recv.size());
        return all_maxloc(send.data(), recv.data(), send.size());
    }

    /*
     * MINLOC CONVENIENCE FUNCTIONS
     */
    template<typename T>
    static bool all_minloc(const T *send, std::pair<T, size_t>* recv, size_t ndata) {
        std::vector<std::pair<T, size_t>> tmp;
        tmp.reserve(ndata);
        for (size_t idata=0ul; idata<ndata; ++idata) tmp.push_back({send[idata], irank()});
        return all_reduce(tmp.data(), recv, MpiMinLoc, ndata);
    }

    template<typename T>
    static bool all_minloc(const T &send, std::pair<T, size_t> &recv) {
        return all_minloc(&send, &recv, 1);
    }

    template<typename T>
    static bool all_minloc(const std::vector<T> &send, std::vector<std::pair<T, size_t>> &recv) {
        ASSERT(send.size()==recv.size());
        return all_minloc(send.data(), recv.data(), send.size());
    }

    template<typename T>
    static bool all_to_all(const T *send, size_t nsend, T *recv, size_t nrecv) {
#ifdef ENABLE_MPI
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
#ifdef ENABLE_MPI
        return MPI_Alltoallv(
                (void *) send, sendcounts, senddispls, mpi_type<T>(),
                (void *) recv, recvcounts, recvdispls, mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        return all_to_all(send, senddispls[0]+sendcounts[0], recv, recvdispls[0]+recvcounts[0]);
#endif
    }

    template<typename T>
    static bool gather(
            const T *send, size_t sendcount, T *recv, size_t recvcount, const size_t& iroot) {
#ifdef ENABLE_MPI
        return MPI_Gather((void *) send, snrw(sendcount), mpi_type<T>(), recv,
                          snrw(recvcount), mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        return all_to_all(send, sendcount, recv, recvcount);
#endif
    }
    template<typename T>
    static bool all_gather(
            const T *send, size_t sendcount, T *recv, size_t recvcount) {
#ifdef ENABLE_MPI
        return MPI_Allgather(
                (void *) send, snrw(sendcount), mpi_type<T>(), recv, snrw(recvcount), mpi_type<T>(), MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        return all_to_all(send, sendcount, recv, recvcount);
#endif
    }



    template<typename T>
    static bool gatherv(
            const T *send, const defs::mpi_count& sendcount,
            T *recv, const defs::mpi_count *recvcounts, const defs::mpi_count *displs, const size_t& iroot) {
#ifdef ENABLE_MPI
        return MPI_Gatherv(
                (void *) send, sendcount, mpi_type<T>(),
                (void *) recv, recvcounts, displs, mpi_type<T>(), iroot, MPI_COMM_WORLD) == MPI_SUCCESS;
#else
        return all_to_all(send, sendcount, recv, recvcounts[0]);
#endif
    }

    template<typename T>
    static bool all_gatherv(
            const T *send, const defs::mpi_count& sendcount,
            T *recv, const defs::mpi_count *recvcounts, const defs::mpi_count *displs) {
#ifdef ENABLE_MPI
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
    static bool gatherv(const T *send, size_t sendcount, T *recv,
                        const defs::inds &recvcounts, const defs::inds &recvdispls, const size_t& iroot) {
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
    static bool gather(const T &send, std::vector<T> &recv, const size_t& iroot) {
        if (recv.size() < nrank()) recv.resize(nrank());
        return gather(&send, 1ul, recv.data(), 1ul, iroot);
    }
    template<typename T>
    static bool all_gather(const T &send, std::vector<T> &recv) {
        if (recv.size() < nrank()) recv.resize(nrank());
        return all_gather(&send, 1ul, recv.data(), 1ul);
    }



    template<typename T>
    static bool gatherv(const std::vector<T> &send, size_t sendcount, std::vector<T> &recv,
                            defs::inds& recvcounts, defs::inds& recvdispls, const size_t& iroot) {
        recvcounts.resize(nrank());
        recvdispls.resize(nrank());
        gather(sendcount, recvcounts, iroot);
        counts_to_displs_consec(recvcounts, recvdispls);
        const size_t nrecv = recvdispls.back() + recvcounts.back();
        if (recv.size() < nrecv) recv.resize(nrecv);
        return gatherv(send.data(), sendcount, recv.data(), recvcounts, recvdispls, iroot);
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
    static bool gatherv(const std::vector<T> &send, std::vector<T> &recv,
                            defs::inds& recvcounts, defs::inds& recvdispls, const size_t& iroot) {
        return gatherv(send, send.size(), recv, recvcounts, recvdispls, iroot);
    }
    template<typename T>
    static bool all_gatherv(const std::vector<T> &send, std::vector<T> &recv,
                            defs::inds& recvcounts, defs::inds& recvdispls) {
        return all_gatherv(send, send.size(), recv, recvcounts, recvdispls);
    }



    template<typename T>
    static bool gatherv(const std::vector<T> &send, size_t sendcount, std::vector<T> &recv, const size_t& iroot) {
        defs::inds recvcounts, recvdispls;
        return gatherv(send, sendcount, recv, recvcounts, recvdispls, iroot);
    }
    template<typename T>
    static bool all_gatherv(const std::vector<T> &send, size_t sendcount, std::vector<T> &recv) {
        defs::inds recvcounts, recvdispls;
        return all_gatherv(send, sendcount, recv, recvcounts, recvdispls);
    }

    template<typename T>
    static bool gatherv(const std::vector<T> &send, std::vector<T> &recv, const size_t& iroot) {
        return gatherv(send, send.size(), recv, iroot);
    }
    template<typename T>
    static bool all_gatherv(const std::vector<T> &send, std::vector<T> &recv) {
        return all_gatherv(send, send.size(), recv);
    }

    static bool i_am(const size_t &i);

    static bool on_node_i_am(const size_t &i);

    static bool i_am_root();

    static bool on_node_i_am_root();

    static bool initialized();

    static bool finalized();

    static void initialize(int *argc = NULL, char ***argv = NULL);

    static void finalize();

    /**
     * call on any rank to terminate the entire communicator
     * @param message
     *  string to output to error log
     */
    static void abort_(std::string message);

    /**
     * call on all ranks to terminate the entire communicator
     * @param message
     *  string to output to error log
     */
    static void abort(std::string message);

};

#endif //M7_MPIWRAPPER_H
