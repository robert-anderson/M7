//
// Created by rja on 23/06/22.
//

#ifndef M7_SMARTPTR_H
#define M7_SMARTPTR_H

#include <memory>

/**
 * since C++11 is the chosen standard for M7, this namespace is included to implement the make_unique and make_shared
 * smart pointer factories, but also includes the ability to create polymorphic instances, which so far is absent even
 * from more recent versions of the standard library.
 *
 * N.B. the "new" C++ keyword should never appear outside of this namespace, and "delete" should appear nowhere in the
 * entire program. Smart pointers were created to automatically manage object lifetimes, there is no good reason to
 * allocate and free objects manually in C++11 and newer.
 */
namespace smart_ptr {
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
     *  "polymorhic" unique pointer
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

#endif //M7_SMARTPTR_H
