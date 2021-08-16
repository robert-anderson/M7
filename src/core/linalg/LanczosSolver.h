//
// Created by rja on 22/10/2020.
//

#ifndef M7_LANCZOSSOLVER_H

#include <src/defs.h>
#include "src/core/linalg/Sparse.h"
#include "src/core/io/Logging.h"

#include <vector>
#include <tuple>
#include <functional>
#include <cassert>
#include <limits>
#include <cmath>
#include <numeric>
#include <random>
#include <algorithm>


namespace lanczos {

    namespace util {

/**
 * @brief Template class to implement positive-definite product
 *
 * ConjugateProduct::prod(a,b) is a function which returns a*b by default.
 * However if the arguments are complex numbers, it returns hermconj(a)*b instead.
 * This structure is required because partial specialization of
 * function template is not allowed in C++.
 */
        template<typename T>
        struct ConjugateProduct {
        public:
            /**
             * @brief Returns a*b for real-type arguments, hermconj(a)*b for complex-type arguments.
             */
            static T prod(T a, T b) { return a * b; }
        };

        template<typename T>
        struct ConjugateProduct<std::complex<T>> {
        public:
            static std::complex<T> prod(std::complex<T> a, std::complex<T> b) {
                return std::conj(a) * b;
            }
        };

/**
 * @brief Returns "mathematical" inner product of v1 and v2.
 *
 * This function is needed because
 * std::inner_product calculates transpose(v1)*v2 instead of dagger(v1)*v2 for complex type.
 *
 */
        template<typename T>
        inline T inner_prod(const std::vector<T> &v1, const std::vector<T> &v2) {
            assert(v1.size() == v2.size());
            return std::inner_product(std::begin(v1), std::end(v1),
                                      std::begin(v2), T(),
                                      [](T a, T b) -> T { return a + b; },
                                      ConjugateProduct<T>::prod);
        }

/**
 * @brief Returns Euclidean norm of given vector.
 */
        template<typename T>
        inline consts::real_t<T> norm(const std::vector<T> &vec) {
            return std::sqrt(std::real(inner_prod(vec, vec)));
            // The norm of any complex vector <v|v> is real by definition.
        }

/**
 * @brief Multiplies each element of vec by a.
 */
        template<typename T1, typename T2>
        inline void scalar_mul(T1 a, std::vector<T2> &vec) {
            for (auto &elem : vec) {
                elem *= a;
            }
        }

/**
 * @brief Normalizes given vector.
 */
        template<typename T>
        inline void normalize(std::vector<T> &vec) {
            scalar_mul(1.0 / norm(vec), vec);
        }


/**
 * @brief Returns 1-norm of given vector.
 */
        template<typename T>
        inline consts::real_t<T> l1_norm(const std::vector<T> &vec) {
            consts::real_t<T> norm = consts::real_t<T>(); // Zero initialization

            for (const T &element : vec) {
                norm += std::abs(element);
            }

            return norm;
        }


/**
 * @brief Orthogonalizes vector `uorth` with respect to orthonormal vectors in `u`.
 *
 * Vectors in `u` must be normalized, but uorth doesn't have to be.
 */
        template<typename T>
        inline void schmidt_orth(std::vector<T> &uorth, const std::vector<std::vector<T>> &u) {
            const auto n = uorth.size();

            for (const auto &uk : u) {
                T innprod = util::inner_prod(uk, uorth);

                for (size_t i = 0; i < n; i++) {
                    uorth[i] -= innprod * uk[i];
                }
            }
        }


/**
 * @brief Returns the significant decimal digits of type T.
 *
 */
        template<typename T>
        inline constexpr int sig_decimal_digit() {
            return (int) (std::numeric_limits<T>::digits *
                          log10(std::numeric_limits<T>::radix));
        }


        template<typename T>
        inline constexpr T minimum_effective_decimal() {
            return pow(10, -sig_decimal_digit<T>());
        }

    }


/**
 * @brief Template class to implement random vector initializer.
 *
 * "Partially specialization of function" is not allowed,
 * so here it is mimicked by wrapping the "init" function with a class template.
 */
    template<typename T>
    struct VectorRandomInitializer {
    public:
        /**
         * @brief Initialize given vector randomly in the range of [-1, 1].
         *
         * For complex type, the real and imaginary part of each element will be initialized in
         * the range of [-1, 1].
         */
        static void init(std::vector<T> &v) {
            std::random_device dev;
            std::mt19937 mt(dev());
            std::uniform_real_distribution<T> rand((T) (-1.0), (T) (1.0));

            size_t n = v.size();
            for (size_t i = 0; i < n; i++) {
                v[i] = rand(mt);
            }
        }
    };


    template<typename T>
    struct VectorRandomInitializer<std::complex<T>> {
    public:
        static void init(std::vector<std::complex<T>> &v) {
            std::random_device dev;
            std::mt19937 mt(dev());
            std::uniform_real_distribution<T> rand((T) (-1.0), (T) (1.0));

            size_t n = v.size();
            for (size_t i = 0; i < n; i++) {
                v[i] = std::complex<T>(rand(mt), rand(mt));
            }
        }
    };


/**
 * @brief Calculation engine for Lanczos algorithm
 */
    template<typename T>
    struct LambdaLanczos {
        /**
         * @brief Matrix-vector multiplication routine.
         *
         * This must be a function to calculate `A*in` and store the result
         * into `out`, where `A` is the matrix to be diagonalized.
         * You can assume the output vector `out` has been initialized with zeros before the `mv_mul` is called.
         */
        std::function<void(const std::vector<T> &in, std::vector<T> &out)> m_mv_mul;

        /** @brief Function to initialize the initial Lanczos vector.
         *
         * After this function called, the output vector will be normalized automatically.
         * Default value is #lambda_lanczos::VectorRandomInitializer::init.
         */
        std::function<void(std::vector<T> &vec)> m_init_vector = VectorRandomInitializer<T>::init;

        /** @brief Dimension of the matrix to be diagonalized. */
        size_t m_matrix_size;
        /** @brief Iteration limit of Lanczos algorithm, set to `matrix_size` automatically. */
        size_t m_max_iteration;
        /** @brief Convergence threshold of Lanczos iteration.
         *
         * `eps` = 1e-12 means that the eigenvalue will be calculated with 12 digits of precision.
         *
         * Default value is system-dependent. On usual 64-bit systems:
         * | type (including complex one)       | size (system-dep.) | `eps`   |
         * | ---------------------------------- | ------------------ | ------- |
         * | float                              | 4 bytes            | 1e-4    |
         * | double                             | 8 bytes            | 1e-12   |
         * | long double                        | 16 bytes           | 1e-19   |
         */
        consts::real_t<T> m_eps = util::minimum_effective_decimal<consts::real_t<T>>() * 1e3;

        /** @brief true to calculate maximum eigenvalue, false to calculate minimum one.*/
        bool m_find_maximum;

        /**
         * @brief Shifts the eigenvalues of the given matrix A.
         *
         * The algorithm will calculate the eigenvalue of matrix (A+`eigenvalue_offset`*E),
         * here E is the identity matrix. The result eigenvalue from `run()` will take this shifting into account,
         * so you don't have to "reshift" the result with `eigenvalue_offset`.
         **/
        consts::real_t<T> m_eigenvalue_offset = 0.0;

        /** @brief (Not necessary to change)
         *
         * Description for those who know the Lanczos algorithm:
         * This is the ratio between the convergence threshold of resulted eigenvalue and the that of
         * tridiagonal eigenvalue. To convergent whole Lanczos algorithm,
         * the convergence threshold for the tridiagonal eigenvalue should be smaller than `eps`.
         */
        consts::real_t<T> m_tridiag_eps_ratio = 1e-1;

        /** @brief (Not necessary to change)
         *
         * This variable specifies the initial reserved size of Lanczos vectors.
         * Controlling this variable might affect reserving efficiency,
         * but that would be very small compared to matrix-vector-multiplication cost.
         */
        size_t m_initial_vector_size = 200;

        /**
          * @brief Constructs Lanczos calculation engine.
          *
          * @param mv_mul Matrix-vector multiplication routine. See #mv_mul for details.
          * @param matrix_size The size of your matrix, i.e. if your matrix is n by n,
          * `matrix_size` should be n.
          * @param find_maximum specifies which of the minimum or maximum eigenvalue to be calculated.
          * By default, `find_maximum=false` so the library will calculates the minimum one.
          */
        LambdaLanczos(std::function<void(const std::vector<T> &, std::vector<T> &)> mv_mul, size_t matrix_size,
                      bool find_maximum = false) :
                m_mv_mul(mv_mul), m_matrix_size(matrix_size), m_max_iteration(matrix_size),
                m_find_maximum(find_maximum) {}

        /**
         * @brief Executes Lanczos algorithm and stores the result into reference variables passed as arguments.
         * @return Lanczos-iteration count
         */
        int run(std::vector<consts::real_t<T>> &eigvalues, std::vector<std::vector<T>> &eigvecs) const {
            const size_t nroot = eigvecs.size();
            assert(eigvalues.size() == nroot);
            assert(0 < m_tridiag_eps_ratio && m_tridiag_eps_ratio < 1);

            std::vector<std::vector<T>> u; // Lanczos vectors
            std::vector<consts::real_t<T>> alpha; // Diagonal elements of an approximated tridiagonal matrix
            std::vector<consts::real_t<T>> beta;  // Subdiagonal elements of an approximated tridiagonal matrix

            const auto n = m_matrix_size;

            u.reserve(m_initial_vector_size);
            alpha.reserve(m_initial_vector_size);
            beta.reserve(m_initial_vector_size);

            u.emplace_back(n, 0.0); // Same as u.push_back(vector<T>(n, 0.0))

            std::vector<T> vk(n, 0.0);

            consts::real_t<T> alphak = 0.0;
            alpha.push_back(alphak);
            consts::real_t<T> betak = 0.0;
            beta.push_back(betak);

            std::vector<T> uk(n);
            m_init_vector(uk);
            util::normalize(uk);
            u.push_back(uk);

            std::vector<consts::real_t<T>> evs(nroot,
                                               0.0); // Calculated eigenvalue (the initial value will never be used)
            std::vector<consts::real_t<T>> pevs(nroot,
                                                std::numeric_limits<consts::real_t<T>>::max()); // Previous eigenvalue

            int itern = m_max_iteration;
            for (size_t k = 1; k <= m_max_iteration; k++) {
                /* vk = (A + offset*E)uk, here E is the identity matrix */
                std::fill(vk.begin(), vk.end(), 0.0);
                m_mv_mul(uk, vk);
                for (size_t i = 0; i < n; i++) {
                    vk[i] += uk[i] * m_eigenvalue_offset;
                }

                alphak = std::real(util::inner_prod(u.back(), vk));

                /* The inner product <uk|vk> is real.
                 * Proof:
                 *     <uk|vk> = <uk|A|uk>
                 *   On the other hand its complex conjugate is
                 *     <uk|vk>^* = <vk|uk> = <uk|A^*|uk> = <uk|A|uk>
                 *   here the condition that matrix A is a symmetric (Hermitian) is used.
                 *   Therefore
                 *     <uk|vk> = <vk|uk>^*
                 *   <uk|vk> is real.
                 */

                alpha.push_back(alphak);

                for (size_t i = 0; i < n; i++) {
                    uk[i] = vk[i] - betak * u[k - 1][i] - alphak * u[k][i];
                }

                util::schmidt_orth(uk, u);

                betak = util::norm(uk);
                beta.push_back(betak);

                for (size_t iroot = 0ul; iroot < nroot; ++iroot)
                    evs[iroot] = find_mth_eigenvalue(alpha, beta, m_find_maximum ? alpha.size() - 2 - iroot : iroot);
                // The first element of alpha is a dummy. Thus its size is alpha.size()-1

                const consts::real_t<T> zero_threshold = util::minimum_effective_decimal<consts::real_t<T>>() * 1e-1;
                if (betak < zero_threshold) {
                    u.push_back(uk);
                    /* This element will never be accessed,
                     * but this "push" guarantees u to always have one more element than
                     * alpha and beta do.
                     */
                    itern = k;
                    break;
                }

                util::normalize(uk);
                u.push_back(uk);

                /*
                 * only break loop if convergence condition is met for all roots
                 */
                bool break_cond = true;
                for (size_t iroot = 0ul; iroot < nroot; ++iroot) {
                    const auto &ev = evs[iroot];
                    const auto &pev = pevs[iroot];
                    if (std::abs(ev - pev) >= std::min(std::abs(ev), std::abs(pev)) * m_eps) {
                        break_cond = false;
                        break;
                    }
                }
                if (break_cond) break;
                else pevs = evs;
            }

            eigvalues = evs;
            for (auto &ev: eigvalues) ev -= m_eigenvalue_offset;

            auto m = alpha.size();
            std::vector<T> cv(m + 1);
            cv[0] = 0.0;
            cv[m] = 0.0;
            cv[m - 1] = 1.0;

            beta[m - 1] = 0.0;

            for (size_t iroot = 0ul; iroot < nroot; ++iroot) {
                auto &eigvec = eigvecs[iroot];
                auto &ev = evs[iroot];

                if (eigvec.size() < n) {
                    eigvec.resize(n);
                }

                for (size_t i = 0; i < n; i++) {
                    eigvec[i] = cv[m - 1] * u[m - 1][i];
                }

                for (size_t k = m - 2; k >= 1; k--) {
                    cv[k] = ((ev - alpha[k + 1]) * cv[k + 1] - beta[k + 1] * cv[k + 2]) / beta[k];

                    for (size_t i = 0; i < n; i++) {
                        eigvec[i] += cv[k] * u[k][i];
                    }
                }
                util::normalize(eigvec);
            }

            return itern;
        }

        int run(consts::real_t<T> &eigvalue, std::vector<T> &eigvec) const {
            std::vector<consts::real_t<T>> eigvalues(1);
            std::vector<std::vector<T>> eigvecs(1);
            auto retval = run(eigvalues, eigvecs);
            eigvalue = eigvalues[0];
            eigvec = std::move(eigvecs[0]);
            return retval;
        }

    private:
        /**
         * @brief Finds the `m`th smaller eigenvalue of given tridiagonal matrix.
         */
        consts::real_t<T> find_mth_eigenvalue(const std::vector<consts::real_t<T>> &alpha,
                                              const std::vector<consts::real_t<T>> &beta,
                                              const size_t m) const {
            consts::real_t<T> eps = m_eps * m_tridiag_eps_ratio;
            consts::real_t<T> mid;
            consts::real_t<T> pmid = std::numeric_limits<consts::real_t<T>>::max();
            consts::real_t<T> r = tridiagonal_eigen_limit(alpha, beta);
            consts::real_t<T> lower = -r;
            consts::real_t<T> upper = r;

            while (upper - lower > std::min(std::abs(lower), std::abs(upper)) * eps) {
                mid = (lower + upper) / 2.0;

                if (num_of_eigs_smaller_than(mid, alpha, beta) >= m + 1) {
                    upper = mid;
                } else {
                    lower = mid;
                }

                if (mid == pmid) {
                    /* This avoids an infinite loop due to zero matrix */
                    break;
                }
                pmid = mid;
            }

            return lower; // The "lower" almost equals the "upper" here.
        }


        /**
         * @brief Computes the upper bound of the absolute value of eigenvalues by Gerschgorin theorem.
         *
         * This routine gives a rough upper bound,
         * but it is sufficient because the bisection routine using
         * the upper bound converges exponentially.
         */
        consts::real_t<T> tridiagonal_eigen_limit(const std::vector<consts::real_t<T>> &alpha,
                                                        const std::vector<consts::real_t<T>> &beta) const {
            consts::real_t<T> r = util::l1_norm(alpha);
            r += 2 * util::l1_norm(beta);

            return r;
        }


        /**
         * @brief Finds the number of eigenvalues of given tridiagonal matrix smaller than `c`.
         *
         * Algorithm from
         * Peter Arbenz et al. / "High Performance Algorithms for Structured Matrix Problems" /
         * Nova Science Publishers, Inc.
         */
        size_t num_of_eigs_smaller_than(consts::real_t<T> c,
                                        const std::vector<consts::real_t<T>> &alpha,
                                        const std::vector<consts::real_t<T>> &beta) const {
            consts::real_t<T> q_i = 1.0;
            size_t count = 0;
            size_t m = alpha.size();

            for (size_t i = 1; i < m; i++) {
                q_i = alpha[i] - c - beta[i - 1] * beta[i - 1] / q_i;
                if (q_i < 0) {
                    count++;
                }
                if (q_i == 0) {
                    q_i = util::minimum_effective_decimal<consts::real_t<T>>();
                }
            }

            return count;
        }
    };


} /* namespace lambda_lanczos */











/**
 * wrapper for the external lambda-lanczos header-only library
 */
struct LanczosSolver {
    std::vector<defs::ham_comp_t> m_evals;
    std::vector<std::vector<defs::wf_t>> m_evecs;

    LanczosSolver(size_t nroot) : m_evals(nroot), m_evecs(nroot) {}

    void solve(const sparse::Matrix<defs::ham_t> &sparse_mat, size_t niter_max) {
        auto mv_mul_fn = [&](const std::vector<defs::wf_t> &in, std::vector<defs::wf_t> &out) {
            sparse_mat.multiply(in, out);
        };
        lanczos::LambdaLanczos<defs::wf_t> external_solver(mv_mul_fn, sparse_mat.nrow());
        external_solver.m_max_iteration = niter_max;
        auto niter = external_solver.run(m_evals, m_evecs);
        log::info("Lanczos found {} roots in {} iterations", m_evals.size(), niter);
    }
};

#endif //M7_LANCZOSSOLVER_H
