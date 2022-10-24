//
// Created by rja on 15/10/22.
//

#ifndef M7_UTIL_POINTER_H
#define M7_UTIL_POINTER_H

#include <iterator>
#include <memory>

namespace ptr {
    template<typename T>
    bool before_begin(const T* ptr, const T* begin) {
        return std::distance(ptr, begin) > 0l;
    }

    template<typename T>
    bool after_begin(const T* ptr, const T* begin) {
        return std::distance(ptr, begin) < 0l;
    }

    template<typename T>
    bool before_end(const T* ptr, const T* end) {
        return std::distance(end, ptr) < 0l;
    }

    template<typename T>
    bool after_end(const T* ptr, const T* end) {
        return std::distance(end, ptr) > 0l;
    }

    /**
     * for a memory block defined by a "begin" pointer and a number of elements n. a pointer is dereferencable to a
     * value within that block if it is equal to or posivitively offset from begin, and is strictly negatively offset
     * from the end pointer: begin + n
     * @tparam T
     *  type of value stored in the block
     * @param ptr
     *  pointer to be checked
     * @param begin
     *  pointer to beginning of the block of memory
     * @param end
     *  pointer to the memory just after the block
     * @return
     *  true if ptr is in the range [begin, end)
     */
    template<typename T>
    bool in_range(const T* ptr, const T* begin, const T* end) {
        return before_end(ptr, end) && !before_begin(ptr, begin);
    }

    template<typename T>
    bool in_range(const T* ptr, const T* begin, uint_t n) {
        return before_end(ptr, begin + n) && !before_begin(ptr, begin);
    }

    /**
     * since C++11 is the chosen standard for M7, this namespace is included to implement the make_unique and
     * make_shared smart pointer factories, but also includes the ability to create polymorphic instances, which so far
     * is absent even from more recent versions of the standard library.
     *
     * N.B. the "new" C++ keyword should never appear outside of this namespace, and "delete" should appear nowhere in
     * the entire program. Smart pointers were created to automatically manage object lifetimes, there is no good reason
     * to allocate and free objects manually in C++11 and newer.
     */
    namespace smart {
        /**
         * @tparam T
         *  both the pointed-to type and the allocated type
         * @tparam Args
         *  types of args to be forwarded to a suitable ctor of T
         * @param args
         *  parameter pack of args to be forwarded to a suitable ctor of T
         * @return
         *  "non-polymorphic" unique pointer
         */
        template<typename T, typename... Args>
        static std::unique_ptr<T> make_unique(Args&&... args) {
            return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
        }
        /**
         * @tparam base_t
         *  pointed-to type
         * @tparam alloc_t
         *  allocated type - must be a subclass of base_t
         * @tparam Args
         *  types of args to be forwarded to a suitable ctor of T
         * @param args
         *  parameter pack of args to be forwarded to a suitable ctor of T
         * @return
         *  "polymorphic" unique pointer
         */
        template<typename base_t, typename alloc_t, typename... Args>
        static std::unique_ptr<base_t> make_poly_unique(Args&&... args) {
            static_assert(std::is_base_of<base_t, alloc_t>::value,
                          "cannot store pointer to allocated object as given base pointer type");
            return std::unique_ptr<base_t>(new alloc_t(std::forward<Args>(args)...));
        }
        /**
         * @tparam T
         *  both the pointed-to type and the allocated type
         * @tparam Args
         *  types of args to be forwarded to a suitable ctor of T
         * @param args
         *  parameter pack of args to be forwarded to a suitable ctor of T
         * @return
         *  "non-polymorphic" shared pointer
         */
        template<typename T, typename... Args>
        static std::shared_ptr<T> make_shared(Args&&... args) {
            return std::shared_ptr<T>(new T(std::forward<Args>(args)...));
        }
        /**
         * @tparam base_t
         *  pointed-to type
         * @tparam alloc_t
         *  allocated type - must be a subclass of base_t
         * @tparam Args
         *  types of args to be forwarded to a suitable ctor of T
         * @param args
         *  parameter pack of args to be forwarded to a suitable ctor of T
         * @return
         *  "polymorhic" shared pointer
         */
        template<typename base_t, typename alloc_t, typename... Args>
        static std::shared_ptr<base_t> make_poly_shared(Args&&... args) {
            static_assert(std::is_base_of<base_t, alloc_t>::value,
                          "cannot store pointer to allocated object as given base pointer type");
            return std::shared_ptr<base_t>(new alloc_t(std::forward<Args>(args)...));
        }
    }
}

#endif //M7_UTIL_POINTER_H
