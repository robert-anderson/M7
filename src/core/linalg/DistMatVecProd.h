//
// Created by rja on 06/01/2022.
//

#ifndef M7_DISTMATVECPROD_H
#define M7_DISTMATVECPROD_H

#include "src/core/parallel/MPIWrapper.h"

template<typename T>
struct DistMatVecProd {
    const size_t m_nelement_local;
    const size_t m_nelement;
    const defs::mpi_counts m_counts;
    const defs::mpi_counts m_displs;
    std::vector<T> m_work_vec;

private:
    defs::mpi_counts make_counts() const {
        defs::mpi_counts counts(mpi::nrank());
        defs::mpi_count tmp = m_nelement_local;
        return mpi::all_gathered(tmp);
    }

public:
    DistMatVecProd(size_t nelement_local):
            m_nelement_local(nelement_local), m_nelement(mpi::all_sum(nelement_local)),
            m_counts(make_counts()), m_displs(mpi::counts_to_displs_consec(m_counts)),
            m_work_vec(nelement_local, 0){}

    void scatter(const T* in, size_t iroot=0ul){
        mpi::scatterv(in, m_nelement, m_work_vec.data(), m_counts.data(), m_displs.data(), iroot);
    }

    void gather(T* out, size_t iroot=0ul){
        mpi::gatherv(m_work_vec.data(), m_nelement, out, m_counts.data(), m_displs.data(), iroot);
    }

    void all_gather(T* out){
        mpi::all_gatherv(m_work_vec.data(), m_nelement, out, m_counts.data(), m_displs.data());
    }

    virtual void scatter_multiply_gather(const T* in, T* out) {

    }
protected:
    virtual void multiply(const T* in, T* out) = 0;
};


#endif //M7_DISTMATVECPROD_H
