//
// Created by rja on 15/06/22.
//

#ifndef M7_TUPLE_H
#define M7_TUPLE_H

#include "utils.h"

namespace utils {
    /**
     * functions acting on Tuples: structures containing diversely-typed data identified only by a constexpr index
     */
    namespace tuple {
        /**
         * execute an operation on each element of a non-const tuple
         * @tparam ind
         *  item index
         * @tparam fn_t
         *  functor type defining a templated call operator which can accept all item_ts as its sole argument
         * @tparam item_ts
         *  parameter pack of all types stored in the tuple
         * @return
         *  void
         */
        template<std::size_t ind = 0, typename fn_t, typename... item_ts>
        inline typename std::enable_if<ind == sizeof...(item_ts), void>::type
        foreach(std::tuple<item_ts...> &, fn_t &) // Unused arguments are given no names.
        {}

        template<std::size_t ind = 0, typename fn_t, typename... item_ts>
        inline typename std::enable_if<ind < sizeof...(item_ts), void>::type
        foreach(std::tuple<item_ts...> &t, fn_t &f) {
            f(std::get<ind>(t));
            foreach<ind + 1, fn_t, item_ts...>(t, f);
        }

        /**
         * execute an operation on each element of a const tuple
         * @tparam ind
         *  item index
         * @tparam fn_t
         *  functor type defining a templated call operator which can accept all item_ts as its sole argument
         * @tparam item_ts
         *  parameter pack of all types stored in the tuple
         * @return
         *  void
         */
        template<std::size_t ind = 0, typename fn_t, typename... item_ts>
        inline typename std::enable_if<ind == sizeof...(item_ts), void>::type
        foreach(const std::tuple<item_ts...> &, fn_t &) // Unused arguments are given no names.
        {}

        template<std::size_t ind = 0, typename fn_t, typename... item_ts>
        inline typename std::enable_if<ind < sizeof...(item_ts), void>::type
        foreach(const std::tuple<item_ts...> &t, fn_t &f) {
            f(std::get<ind>(t));
            foreach<ind + 1, fn_t, item_ts...>(t, f);
        }


        /*
         * modifiable / modifiable
         */
        /**
         * @tparam ind
         * @tparam fn_t
         * @tparam item_ts
         * @param t1
         * @param t2
         * @param f
         * @return
         */
        template<std::size_t ind = 0, typename fn_t, typename... item_ts>
        inline typename std::enable_if<ind == sizeof...(item_ts), void>::type
        foreach(std::tuple<item_ts...> &t1, std::tuple<item_ts...> &t2, fn_t &f){}

        template<std::size_t ind = 0, typename fn_t, typename... item_ts>
        inline typename std::enable_if<ind < sizeof...(item_ts), void>::type
        foreach(std::tuple<item_ts...> &t1, std::tuple<item_ts...> &t2, fn_t &f) {
            f(std::get<ind>(t1), std::get<ind>(t2));
            foreach<ind + 1, fn_t, item_ts...>(t1, t2, f);
        }


        /*
         * modifiable / const
         */
        template<std::size_t ind = 0, typename fn_t, typename... item_ts>
        inline typename std::enable_if<ind == sizeof...(item_ts), void>::type
        foreach(std::tuple<item_ts...> &, const std::tuple<item_ts...> &, fn_t &) {}

        template<std::size_t ind = 0, typename fn_t, typename... item_ts>
        inline typename std::enable_if<ind < sizeof...(item_ts), void>::type
        foreach(std::tuple<item_ts...> &t1, const std::tuple<item_ts...> &t2, fn_t &f) {
            f(std::get<ind>(t1), std::get<ind>(t2));
            foreach<ind + 1, fn_t, item_ts...>(t1, t2, f);
        }

        /*
         * const / const
         */
        template<std::size_t ind = 0, typename fn_t, typename... item_ts>
        inline typename std::enable_if<ind == sizeof...(item_ts), void>::type
        foreach(const std::tuple<item_ts...> &, const std::tuple<item_ts...> &, fn_t &) {}

        template<std::size_t ind = 0, typename fn_t, typename... item_ts>
        inline typename std::enable_if<ind < sizeof...(item_ts), void>::type
        foreach(const std::tuple<item_ts...> &t1, const std::tuple<item_ts...> &t2, fn_t &f) {
            f(std::get<ind>(t1), std::get<ind>(t2));
            foreach<ind + 1, fn_t, item_ts...>(t1, t2, f);
        }
    }
}

#endif //M7_TUPLE_H
