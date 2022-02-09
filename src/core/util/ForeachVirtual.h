//
// Created by anderson on 2/9/22.
//

#ifndef M7_FOREACHVIRTUAL_H
#define M7_FOREACHVIRTUAL_H

#include "utils.h"

namespace foreach_virtual {

    class ExitLoop: public std::exception {
        virtual const char* what() const throw() {
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
            const size_t m_nterm;
            inds_t<nind> m_inds;
            const inds_t<nind>& inds() const {return m_inds;}
            const size_t& operator[](const size_t& i){
                ASSERT(i<nind);
                return m_inds[i];
            }
            size_t sum() const {
                return std::accumulate(m_inds.cbegin(), m_inds.cend(), 0);
            }
            virtual void body(const inds_t<nind>&) = 0;

            Base(size_t nterm) : m_nterm(nterm) {}

            virtual void loop() = 0;

        };

        template<size_t nind>
        struct Unrestricted : Base<nind> {
            using Base<nind>::m_inds;
            using Base<nind>::body;
            const inds_t<nind> m_shape;
        private:
            static size_t nterm(const inds_t<nind> &shape) {
                if (!nind) return 0;
                size_t n = 1;
                for (size_t i = 0ul; i != nind; ++i) n *= shape[i];
                return n;
            }

        public:
            Unrestricted(const inds_t<nind> &shape) : Base<nind>(nterm(shape)), m_shape(shape) {}
            Unrestricted(const size_t& extent) : Unrestricted(array_utils::filled<size_t, nind>(extent)){}

            template<size_t ilevel>
            void level_loop(tags::Ind<ilevel>) {
                constexpr size_t iind = ilevel - 1;
                auto &ind = Base<nind>::m_inds[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) {
                    try {level_loop(tags::Ind<ilevel + 1>());}
                    catch (const ExitLoop& ex) {throw ex;}
                }
            }

            void level_loop(tags::Ind<nind>) {
                constexpr size_t iind = nind - 1;
                auto &ind = m_inds[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) {
                    try {body(m_inds);}
                    catch (const ExitLoop& ex) {throw ex;}
                }
            }

            void top_loop(tags::Bool<true>) {}

            void top_loop(tags::Bool<false>) {
                try {level_loop(tags::Ind<1>());}
                catch (const ExitLoop& ex) {}

            }

            void loop() override {
                top_loop(tags::Bool<nind == 0>());
            }
        };

        template<size_t nind, bool strict = true, bool ascending = true>
        struct Ordered : Base<nind> {
            using Base<nind>::m_inds;
            using Base<nind>::body;
            using inds_t = std::array<size_t, nind>;
            const size_t m_n;

            static size_t nterm(size_t n) {
                if (!nind) return 0ul;
                return integer_utils::combinatorial(strict ? n : (n + nind) - 1, nind);
            }

        public:
            Ordered(size_t n) :
                    Base<nind>(nterm(n)), m_n(n) {}

            template<size_t ilevel>
            void level_loop(tags::Ind<ilevel>) {
                constexpr size_t iind = ascending ? (nind - ilevel) : (ilevel - 1);
                constexpr size_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_inds[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_inds[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) {
                    try {level_loop(tags::Ind<ilevel + 1>());}
                    catch (const ExitLoop& ex) {throw ex;}
                }
            }

            void level_loop(tags::Ind<nind>) {
                constexpr size_t iind = ascending ? 0 : nind - 1;
                constexpr size_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_inds[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_inds[ascending ? iind + 1 : iind - 1] + !strict;
                for (ind = 0ul; ind < extent; ++ind) {
                    try {body(m_inds);}
                    catch (const ExitLoop& ex) {throw ex;}
                }
            }

            void top_loop(tags::Bool<true>) {}

            void top_loop(tags::Bool<false>) {
                try {level_loop(tags::Ind<1>());}
                catch (const ExitLoop& ex) {}
            }

            void loop() override {
                top_loop(tags::Bool<nind == 0>());
            }
        };
    }

    /**
     * "run time number of dimensions"
     */
    namespace rtnd {
        using inds_t = std::vector<size_t>;
        using body_fn_t = std::function<void()>;
        using cr_body_fn_t = const body_fn_t &;

        struct Base {
            const size_t m_nind;
            const size_t m_nterm;
            inds_t m_inds;

            const inds_t& inds() const {return m_inds;}

            const size_t& operator[](const size_t& i){
                ASSERT(i<m_nind);
                return m_inds[i];
            }

            size_t ind_sum() const {
                return std::accumulate(m_inds.cbegin(), m_inds.cend(), 0);
            }

            virtual void operator()(cr_body_fn_t body) = 0;

            template<typename ...Args>
            void operator()(cr_body_fn_t body, Base& next, Args&... rest) {
                auto this_body = [&](){
                    next(body, rest...);
                };
                (*this)(this_body);
            }

            Base(size_t nind, size_t nterm) : m_nind(nind),
                                              m_nterm(nterm), m_inds(m_nind, 0) {}

        };

        template<typename ...Args>
        static void chain(cr_body_fn_t body, Base& first, Args&... rest) {
            first(body, rest...);
        }

        struct Unrestricted : Base {
            const inds_t m_shape;
            using Base::operator();
        private:
            static size_t nterm(const inds_t &shape) {
                if (shape.empty()) return 0;
                size_t n = 1;
                for (const auto &i: shape) n *= i;
                return n;
            }

        public:
            Unrestricted(const inds_t &shape) : Base(shape.size(), nterm(shape)), m_shape(shape) {}
            Unrestricted(const size_t& nind, const size_t& extent) : Unrestricted(std::vector<size_t>(nind, extent)){}

            void level_loop(cr_body_fn_t body, size_t ilevel) {
                const auto &iind = ilevel - 1;
                auto &ind = m_inds[iind];
                const auto &extent = m_shape[iind];
                if (ilevel < m_nind)
                    for (ind = 0ul; ind < extent; ++ind) level_loop(body, ilevel + 1);
                else
                    for (ind = 0ul; ind < extent; ++ind) body();
            }

            void operator()(cr_body_fn_t body) override {
                if (m_nind == 0) return;
                level_loop(body, 1);
            }
        };

        template<bool strict = true, bool ascending = true>
        struct Ordered : Base {
            const size_t m_n;
            using Base::operator();

            static size_t nterm(size_t n, size_t r) {
                if (!r) return 0ul;
                return integer_utils::combinatorial(strict ? n : (n + r) - 1, r);
            }

        public:
            Ordered(const size_t &n, const size_t &r) :
                    Base(r, nterm(n, r)), m_n(n) {}

            void level_loop(cr_body_fn_t body, size_t ilevel) {
                const size_t iind = ascending ? (m_nind - ilevel) : (ilevel - 1);
                const size_t iind_unrestrict = ascending ? m_nind - 1 : 0;
                auto &ind = m_inds[iind];
                const auto extent = iind == iind_unrestrict ? m_n : m_inds[ascending ? iind + 1 : iind - 1] + !strict;
                if (ilevel < m_nind)
                    for (ind = 0ul; ind < extent; ++ind) level_loop(body, ilevel + 1);
                else
                    for (ind = 0ul; ind < extent; ++ind) body();
            }

            void operator()(cr_body_fn_t body) override {
                if (m_nind == 0) return;
                level_loop(body, 1);
            }
        };
    }
}


#endif //M7_FOREACHVIRTUAL_H
