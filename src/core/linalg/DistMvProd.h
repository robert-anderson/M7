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
                m_in_vec(m_nelement, 0), m_out_vec(nelement_local, 0) {}

        T* get_in_ptr(const T *in) const {
            auto ptr = in;
            if (!ptr) ptr = m_in_vec.data();
            return const_cast<T*>(ptr);
        }

        void bcast(const T *in) {
            mpi::bcast(get_in_ptr(in), m_nelement);
        }

        void gather(T *out) {
            mpi::gatherv(m_out_vec.data(), m_nelement_local, out, m_counts.data(), m_displs.data(), 0ul);
        }

        void all_gather(T *out) {
            mpi::all_gatherv(m_out_vec.data(), m_nelement_local, out, m_counts.data(), m_displs.data());
        }

        void parallel_multiply(const T *in, T *out, bool all_out = false) {
            bcast(in);
            multiply(in);
            all_out ? all_gather(out) : gather(out);
        }

        void parallel_multiply(const std::vector<T> &in, std::vector<T> &out, bool all_out = false) {
            if (all_out || mpi::i_am_root()) out.resize(m_nelement);
            parallel_multiply(in.data(), out.data(), all_out);
        }

    protected:
        virtual void multiply(const T *in) = 0;
    };

    template<typename T>
    struct Sparse: public Base<T> {
        const sparse::Matrix<T>& m_mat;
        Sparse(const sparse::Matrix<T>& mat): Base<T>(mat.nrow()), m_mat(mat){}

    protected:
        using Base<T>::m_in_vec;
        using Base<T>::m_out_vec;
        void multiply(const T *in) override {
            m_mat.multiply(Base<T>::get_in_ptr(in), m_out_vec.data());
        }
    };

    /*
    struct OnTheFly {

    };
     */
}


#endif //M7_DISTMVPROD_H
