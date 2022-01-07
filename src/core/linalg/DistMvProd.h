//
// Created by rja on 06/01/2022.
//

#ifndef M7_DISTMVPROD_H
#define M7_DISTMVPROD_H

#include "src/core/parallel/MPIWrapper.h"
#include "Sparse.h"

namespace dist_mv_prod {
    template<typename T>
    struct Base {
        const size_t m_nelement_local;
        const size_t m_nelement;
        const defs::mpi_counts m_counts;
        const defs::mpi_counts m_displs;
        std::vector<T> m_in_vec;
        std::vector<T> m_out_vec;

    private:
        defs::mpi_counts make_counts() const {
            defs::mpi_counts counts(mpi::nrank());
            defs::mpi_count tmp = m_nelement_local;
            return mpi::all_gathered(tmp);
        }

    public:
        Base(size_t nelement_local) :
                m_nelement_local(nelement_local), m_nelement(mpi::all_sum(nelement_local)),
                m_counts(make_counts()), m_displs(mpi::counts_to_displs_consec(m_counts)),
                m_in_vec(nelement_local, 0), m_out_vec(nelement_local, 0) {}

        void scatter(const T *in, size_t iroot = 0ul) {
            mpi::scatterv(in, m_nelement, m_in_vec.data(), m_counts.data(), m_displs.data(), iroot);
        }

        void gather(T *out, size_t iroot = 0ul) {
            mpi::gatherv(m_out_vec.data(), m_nelement, out, m_counts.data(), m_displs.data(), iroot);
        }

        void all_gather(T *out) {
            mpi::all_gatherv(m_out_vec.data(), m_nelement, out, m_counts.data(), m_displs.data());
        }

        void parallel_multiply(const T *in, T *out, bool all_out = false) {
            scatter(in);
            multiply();
            all_out ? all_gather(out) : gather(out);
        }

    protected:
        virtual void multiply() = 0;
    };

    template<typename T>
    struct Sparse: Base<T> {
        const sparse::Matrix<T>& m_mat;
        Sparse(const sparse::Matrix<T>& mat): Base<T>(mat.nrow()), m_mat(mat){}

    protected:
        void multiply() override {
            m_mat.multiply(m_in_vec, m_out_vec);
        }
    };

    /*
    struct OnTheFly {

    };
     */
}


#endif //M7_DISTMVPROD_H
