//
// Created by Robert J. Anderson on 2/9/22.
//

#ifndef M7_FOREACHVIRTUAL_H
#define M7_FOREACHVIRTUAL_H


#include <utility>

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/util/utils.h>



namespace foreach_virtual {

    class ExitLoop : public std::exception {
        virtual const char *what() const throw() {
            return "Loop body requested early termination of loop";
        }
    };
    /**
     * "compile time number of dimensions"
     */
    namespace ctnd {

        template<size_t nind>
        using inds_t = std::array<size_t, nind>;

        template<size_t nind>
        struct Base {
        protected:
            static constexpr size_t c_nind = nind;
            /**
             * work space into which each set of indices is inserted: should not be directly accessed in derived classes
             */
            inds_t<nind> m_value;
            /**
             * integer index of the current iteration (body call)
             */
            size_t m_iiter;
            /**
             * length of the index iterator
             */
            const size_t m_niter;
        public:
            /**
             * getter for the current index value so that body implementations make immutable reference to the indices
             */
            const inds_t<nind> &value() const { return m_value; }

            /**
             * getter for the current iteration index
             */
            const size_t &iiter() const { return m_iiter; }

            /**
             * getter for the length of the iterator
             */
            const size_t &niter() const { return m_niter; }

            /**
             * function to be called each time a new set of indices is formed
             */
            virtual void body(const inds_t <nind> &value, size_t iiter) = 0;

        protected:
            /**
             * function defining the algorithm by which sets of indices are generated
             */
            virtual void throwing_loop() = 0;

        public:

            Base(size_t niter): m_niter(niter){}

            void loop() {
                if (nind) {
                    try {
                        throwing_loop();
                        DEBUG_ASSERT_EQ(m_iiter+1, m_niter, "loop completed after incorrect number of iterations");
                    }
                    catch (const ExitLoop &) {}
                } else m_iiter = 0ul;
            }
        };

        template<size_t nind>
        struct Unrestricted : Base<nind> {
            using Base<nind>::m_value;
            using Base<nind>::m_iiter;
            using Base<nind>::m_niter;
            using Base<nind>::body;
            const inds_t<nind> m_shape;
        private:
            /**
             * @param shape
             *  extents of each index
             * @return
             *  expected number of calls to body() without early termination
             */
            static size_t nterm(const inds_t<nind> &shape) {
                if (!nind) return 0;
                size_t n = 1;
                for (size_t i = 0ul; i != nind; ++i) n *= shape[i];
                return n;
            }

            /**
             * loop through the values of the ilevel element of m_value between 0 and the extent of that dimension
             * @tparam ilevel
             *  current level (position in the m_value array)
             */
            template<size_t ilevel>
            void level_loop(tags::Int<ilevel>) {
                constexpr size_t iind = ilevel - 1;
                auto &ind = Base<nind>::m_value[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) {
                    try { level_loop(tags::Int<ilevel + 1>()); }
                    catch (const ExitLoop &ex) { throw ex; }
                }
            }

            /**
             * overload for when the last index has been reached
             */
            void level_loop(tags::Int<nind>) {
                constexpr size_t iind = nind - 1;
                auto &ind = m_value[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) {
                    ++m_iiter;
                    try { body(m_value, m_iiter); }
                    catch (const ExitLoop &ex) { throw ex; }
                }
            }

            /**
             * in the edge case that the nind is 0, do nothing but set the iteration counter
             */
            void top_loop(tags::Int<true>) {
                m_iiter = 0ul;
            }

            /**
             * if nind is nonzero, start at the first index
             */
            void top_loop(tags::Int<false>) {
                try { level_loop(tags::Int<1>()); }
                catch (const ExitLoop &ex) {throw ex;}
            }

        public:

            Unrestricted(const inds_t<nind> &shape) : Base<nind>(nterm(shape)), m_shape(shape){}

            Unrestricted(size_t extent=0ul) : Unrestricted(array_utils::filled<size_t, nind>(extent)) {}


        protected:
            void throwing_loop() override {
                m_iiter = ~0ul;
                top_loop(tags::Int<nind == 0>());
            }
        };

        template<size_t nind, bool strict = true, bool ascending = true>
        struct Ordered : Base<nind> {
            using Base<nind>::m_value;
            using Base<nind>::m_iiter;
            using Base<nind>::m_niter;
            using Base<nind>::body;
            using inds_t = std::array<size_t, nind>;
        protected:
            size_t m_n;

        private:
            static size_t nterm(size_t n) {
                if (!nind) return 0ul;
                return integer_utils::combinatorial(strict ? n : (n + nind) - 1, nind);
            }

            template<size_t ilevel>
            void level_loop(tags::Int<ilevel>) {
                constexpr size_t iind = ascending ? (nind - ilevel) : (ilevel - 1);
                constexpr size_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) {
                    try { level_loop(tags::Int<ilevel + 1>()); }
                    catch (const ExitLoop &ex) { throw ex; }
                }
            }

            void level_loop(tags::Int<nind>) {
                constexpr size_t iind = ascending ? 0 : nind - 1;
                constexpr size_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) {
                    ++m_iiter;
                    try { body(m_value, m_iiter); }
                    catch (const ExitLoop &ex) { throw ex; }
                }
            }

            void top_loop(tags::Int<true>) {}

            void top_loop(tags::Int<false>) {
                try { level_loop(tags::Int<1>()); }
                catch (const ExitLoop &ex) {throw ex;}
            }

        public:

            Ordered(size_t n) : Base<nind>(nterm(n)), m_n(n){}

        protected:
            void throwing_loop() override {
                m_iiter = ~0ul;
                top_loop(tags::Int<nind == 0>());
            }
        };

        /**
         * functional objects are convenient but entail non-negligible overhead. If this is non-critical, the following
         * definitions may be confidently used
         */
        namespace lambda {
            template<typename foreach_t>
            struct Lambda : foreach_t {
                typedef ctnd::inds_t<foreach_t::c_nind> inds_t;
                typedef std::function<void(const inds_t &value, size_t iiter)> body_fn_t;
                body_fn_t m_body_fn;
                template<typename ...Args>
                Lambda(body_fn_t fn, Args&&... args): foreach_t(std::forward<Args...>(args...)), m_body_fn(fn){}
                void body(const inds_t &value, size_t iiter) override {
                    m_body_fn(value, iiter);
                }
            };
            template<size_t nind>
            struct Unrestricted : Lambda<ctnd::Unrestricted<nind>> {
                typedef Lambda<ctnd::Unrestricted<nind>> base_t;
                using typename base_t::body_fn_t;
                using typename base_t::inds_t;
                Unrestricted(body_fn_t fn, inds_t shape): base_t(fn, std::move(shape)){}
                Unrestricted(body_fn_t fn, size_t extent): base_t(fn, extent){}
            };
            template<size_t nind, bool strict = true, bool ascending = true>
            using Ordered = Lambda<ctnd::Ordered<nind, strict, ascending>>;
        }
    }

    /**
     * "run time number of dimensions"
     */
    namespace rtnd {
        using inds_t = std::vector<size_t>;

        struct Base {
        protected:
            /**
             * work space into which each set of indices is inserted: should not be directly accessed in derived classes
             */
            inds_t m_value;
            /**
             * integer index of the current iteration (body call)
             */
            size_t m_iiter;
            /**
             * length of the index iterator
             */
            const size_t m_niter;
        public:
            /**
             * number of dimensions: length of the index array
             */
            const size_t m_nind;

            const inds_t &value() const;

            const size_t &iiter() const;

            const size_t &niter() const;

            Base(size_t nind, size_t niter);

            /**
             * function to be called each time a new set of indices is formed
             */
            virtual void body(const inds_t &value, size_t iiter) = 0;

            /**
             * function with defines the looping logic, calls body, and increments the iteration counter
             */
            virtual void throwing_loop() = 0;

            /**
             * exposed wrapper function which catches any instance of the ExitLoop exception thrown, and gracefully
             * handles the zero-dimensional edge case
             */
            void loop();
        };

        struct Unrestricted : Base {
            const inds_t m_shape;
        private:
            static size_t nterm(const inds_t &shape);

            void level_loop(size_t ilevel);
        public:
            explicit Unrestricted(inds_t shape);

            Unrestricted(size_t nind, size_t extent);

            void throwing_loop() override;
        };

        template<bool strict = true, bool ascending = true>
        struct Ordered : Base {
        private:
            const size_t m_n;

            static size_t nterm(size_t n, size_t r) {
                if (!r) return 0ul;
                return integer_utils::combinatorial(strict ? n : (n + r) - 1, r);
            }

            void level_loop(size_t ilevel) {
                const size_t iind = ascending ? (m_nind - ilevel) : (ilevel - 1);
                const size_t iind_unrestrict = ascending ? m_nind - 1 : 0;
                auto &ind = m_value[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_value[ascending ? iind + 1 : iind - 1] + !strict;
                try {
                    if (ilevel < m_nind)
                        for (ind = 0ul; ind < extent; ++ind) level_loop(ilevel + 1);
                    else {
                        for (ind = 0ul; ind < extent; ++ind) {
                            ++m_iiter;
                            body(m_value, m_iiter);
                        }
                    }
                }
                catch (const ExitLoop& ex){throw ex;}
            }

        public:
            Ordered(size_t n, size_t r) : Base(r, nterm(n, r)), m_n(n) {}

            void throwing_loop() override {
                m_iiter = ~0ul;
                level_loop(1);
            }
        };

        /**
         * functional objects are convenient but entail non-negligible overhead. If this is non-critical, the following
         * definitions may be confidently used
         */
        namespace lambda {
            template<typename foreach_t>
            struct Lambda : foreach_t {
                typedef std::function<void(const defs::inds &value, size_t iiter)> body_fn_t;
                body_fn_t m_body_fn;
                template<typename ...Args>
                Lambda(body_fn_t fn, Args&&... args):
                    foreach_t(std::forward<Args>(args)...), m_body_fn(std::move(fn)){}
                void body(const inds_t &value, size_t iiter) override {
                    m_body_fn(value, iiter);
                }
            };
            struct Unrestricted : Lambda<rtnd::Unrestricted> {
                typedef Lambda<rtnd::Unrestricted> base_t;
                Unrestricted(body_fn_t fn, inds_t shape): base_t(std::move(fn), std::move(shape)){}
                Unrestricted(body_fn_t fn, size_t nind, size_t extent): base_t(std::move(fn), nind, extent){}
            };

            template<bool strict = true, bool ascending = true>
            struct Ordered: Lambda<rtnd::Ordered<strict, ascending>> {
                typedef Lambda<rtnd::Ordered<strict, ascending>> base_t;
                using typename base_t::body_fn_t;
                Ordered(body_fn_t fn, size_t n, size_t r): base_t(std::move(fn), n, r){}
            };
        }
    }
}


#endif //M7_FOREACHVIRTUAL_H
