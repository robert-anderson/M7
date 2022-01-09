//
// Created by rja on 06/01/2022.
//

#ifndef M7_DISTMVPROD_H
#define M7_DISTMVPROD_H

#include "src/core/parallel/MPIWrapper.h"
#include "Sparse.h"

namespace dist_mv_prod {
    /**
     * matrix representation-agnostic base class for computing the product Mv over all MPI ranks
     * @tparam T
     *  type of the elements of v and Mv
     */
    template<typename T>
    struct Base {
        /**
         * number of elements of the product vector (Mv) stored on this rank
         * this is the same as the number of rows of the matrix multiplied by the input vector v on this rank
         */
        const size_t m_nrow_local;
        /**
         * total number of elements in the product vector (Mv) over all ranks
         */
        const size_t m_nrow;
        /**
         * m_nrow_local values for all ranks
         */
        const defs::mpi_counts m_counts;
        /**
         * cumulative counts
         */
        const defs::mpi_counts m_displs;
        /**
         * input vector buffer whose size is given by the number of columns of M. This size is not assumed to be known
         * until an input vector is passed to the multiply method, whereupon all m_in_vecs are resized accordingly
         */
        std::vector<T> m_v;
        /**
         * partial product vector. the full Mv is recovered when these are combined consecutively across all ranks
         */
        std::vector<T> m_partial_mv;

    private:
        defs::mpi_counts make_counts() const {
            defs::mpi_counts counts(mpi::nrank());
            defs::mpi_count tmp = m_nrow_local;
            return mpi::all_gathered(tmp);
        }

    public:
        Base(size_t nelement_mv_local) :
                m_nrow_local(nelement_mv_local), m_nrow(mpi::all_sum(nelement_mv_local)),
                m_counts(make_counts()), m_displs(mpi::counts_to_displs_consec(m_counts)),
                m_v(mpi::i_am_root() ? 0 : m_nrow, 0), m_partial_mv(nelement_mv_local, 0) {}
                
        void parallel_multiply(const T *v, size_t v_size, T *mv, bool all_gather_mv = false) {
            bcast(v, v_size);
            multiply(v, v_size);
            all_gather_mv ? all_gather(mv) : gather(mv);
        }

        void parallel_multiply(const std::vector<T> &v, std::vector<T> &mv, bool all_gather_mv = false) {
            if (all_gather_mv || mpi::i_am_root()) mv.resize(m_nrow);
            parallel_multiply(v.data(), v.size(), mv.data(), all_gather_mv);
        }
        
    protected:
        T* get_v_ptr(const T *v, size_t v_size) {
            auto ptr = v;
            if (!ptr) {
                m_v.resize(v_size);
                ptr = m_v.data();
            }
            return const_cast<T*>(ptr);
        }

        void bcast(const T *v, size_t v_size) {
            mpi::bcast(get_v_ptr(v, v_size), m_nrow);
        }

        void gather(T *mv) {
            REQUIRE_TRUE_ALL(mv || !mpi::i_am_root(), "gathering pointer must be non-null");
            mpi::gatherv(m_partial_mv.data(), m_nrow_local, mv, m_counts.data(), m_displs.data(), 0ul);
        }

        void all_gather(T *mv) {
            REQUIRE_TRUE_ALL(mv, "all gathering pointers must be non-null");
            mpi::all_gatherv(m_partial_mv.data(), m_nrow_local, mv, m_counts.data(), m_displs.data());
        }

        virtual void multiply(const T *v, size_t v_size) = 0;
    };

    template<typename T>
    struct Sparse: public Base<T> {
        const sparse::Matrix<T>& m_mat;
        Sparse(const sparse::Matrix<T>& mat): Base<T>(mat.nrow()), m_mat(mat){}

    protected:
        void multiply(const T *v, size_t v_size) override {
            m_mat.multiply(Base<T>::get_v_ptr(v, v_size), Base<T>::m_partial_mv.data());
        }
    };

    /*
    struct OnTheFly {

    };
     */
}


#endif //M7_DISTMVPROD_H
