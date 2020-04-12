//
// Created by Robert John Anderson on 2020-04-12.
//

#ifndef M7_ALIGNEDALLOCATOR_H
#define M7_ALIGNEDALLOCATOR_H

#include <cstddef>
#include <assert.h>
#include <limits>
#include <src/defs.h>

template<class T, size_t alignment>
class AlignedAllocator {
    T *m_malloc_ptr = nullptr;
public:
    // type definitions
    typedef T value_type;
    typedef T *pointer;
    typedef const T *const_pointer;
    typedef T &reference;
    typedef const T &const_reference;
    typedef std::size_t size_type;

    // rebind allocator to type U
    template<class U>
    struct rebind {
        typedef AlignedAllocator<U, alignment> other;
    };

    // return address of values
    pointer address(reference value) const {
        return &value;
    }

    const_pointer address(const_reference value) const {
        return &value;
    }

    /* constructors and destructor
     * - nothing to do because the allocator has no state
     */
    AlignedAllocator() throw() {
        //assert(alignof(T) == alignment);
    }

    AlignedAllocator(const AlignedAllocator &alloc) throw() {
        m_malloc_ptr = alloc.m_malloc_ptr;
    }

    template<class U>
    AlignedAllocator(const AlignedAllocator<U, alignment> &) throw() {
    }

    ~AlignedAllocator() throw() {
    }

    // return maximum number of elements that can be allocated
    size_type max_size() const throw() {
        return std::numeric_limits<std::size_t>::max() / sizeof(T);
    }

    // allocate but don't initialize num elements of type T
    pointer allocate(size_type num, const void * = 0) {
        m_malloc_ptr = (pointer) (::operator new(num * sizeof(T) + alignment));
        const auto overshoot = ((size_t) m_malloc_ptr) % alignment;
        return pointer((char *) m_malloc_ptr + alignment - overshoot);
    }

    // initialize elements of allocated storage p with value value
    void construct(pointer p, const T &value) {
        // initialize memory with placement new
        new((void *) p)T(value);
    }

    // destroy elements of initialized storage p
    void destroy(pointer p) {
        // destroy objects by calling their destructor
        p->~T();
    }

    // deallocate storage p of deleted elements
    void deallocate(pointer p, size_type num) {
        ::operator delete((void *) m_malloc_ptr);
    }
};

/*
// return that all specializations of this allocator are interchangeable
template<class T1, class T2, size_t alignment>
bool operator==(const MyAlloc<T1, alignment> &,
                const MyAlloc<T2, alignment> &) throw() {
    return true;
}

template<class T1, class T2, size_t alignment>
bool operator!=(const MyAlloc<T1, alignment> &,
                const MyAlloc<T2, alignment> &) throw() {
    return false;
}
 */


#endif //M7_ALIGNEDALLOCATOR_H
