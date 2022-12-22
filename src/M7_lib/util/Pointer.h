//
// Created by rja on 15/10/22.
//

#ifndef M7_UTIL_POINTER_H
#define M7_UTIL_POINTER_H

#include <iterator>
#include <memory>

namespace ptr {
    template<typename T>
    const T* before_begin(const T* ptr, const T* begin) {
        return (std::distance(ptr, begin) > 0l) ? ptr : nullptr;
    }
    template<typename T>
    T* before_begin(T* ptr, const T* begin) {
        return const_cast<T*>(before_begin(const_cast<const T*>(ptr), begin));
    }

    template<typename T>
    const T* after_begin(const T* ptr, const T* begin) {
        return (std::distance(ptr, begin) < 0l) ? ptr : nullptr;
    }
    template<typename T>
    T* after_begin(T* ptr, const T* begin) {
        return const_cast<T*>(after_begin(const_cast<const T*>(ptr), begin));
    }

    template<typename T>
    const T* before_end(const T* ptr, const T* end) {
        return (std::distance(end, ptr) < 0l) ? ptr : nullptr;
    }
    template<typename T>
    T* before_end(T* ptr, const T* end) {
        return const_cast<T*>(before_end(const_cast<const T*>(ptr), end));
    }

    template<typename T>
    const T* after_end(const T* ptr, const T* end) {
        return (std::distance(end, ptr) > 0l) ? ptr : nullptr;
    }
    template<typename T>
    T* after_end(T* ptr, const T* end) {
        return const_cast<T*>(after_end(const_cast<const T*>(ptr), end));
    }

    /**
     * for a memory block defined by a "begin" pointer and a number of elements n. a pointer is dereferencable to a
     * value within that block if it is equal to or positively offset from begin, and is strictly negatively offset
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
     *  ptr if ptr is in the range [begin, end), else nullptr
     */
    template<typename T>
    T* in_range(T* ptr, const T* begin, const T* end) {
        return (before_end(ptr, end) && !before_begin(ptr, begin)) ? ptr : nullptr;
    }
    template<typename T>
    const T* in_range(const T* ptr, const T* begin, const T* end) {
        return const_cast<const T*>(in_range(const_cast<T*>(ptr), begin, end));
    }

    template<typename T>
    T* in_range(T* ptr, const T* begin, uint_t n) {
        return in_range(ptr, begin, begin + n);
    }
    template<typename T>
    const T* in_range(const T* ptr, const T* begin, uint_t n) {
        return in_range(ptr, begin, begin + n);
    }

    /**
     * update the pointer at ptr but return its pre-update value
     */
    template<typename T>
    T* pre_update(T*& ptr, T* new_ptr) {
        auto tmp = ptr;
        ptr = new_ptr;
        return tmp;
    }

    template<typename T>
    const T* pre_update(const T*& ptr, const T* new_ptr) {
        return const_cast<const T*>(pre_update(const_cast<T*&>(ptr), const_cast<T*>(new_ptr)));
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
