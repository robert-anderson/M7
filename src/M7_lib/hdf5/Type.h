//
// Created by Robert J. Anderson on 13/12/2020.
//

#ifndef M7_HDF5_TYPE_H
#define M7_HDF5_TYPE_H

#include <string>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <M7_lib/parallel/MPIWrapper.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/nd/NdFormat.h>

namespace hdf5 {
#ifdef H5_HAVE_PARALLEL
    constexpr bool c_have_parallel = true;
#else
    constexpr bool c_have_parallel = false;
#endif

    static_assert(c_have_parallel, "HDF5 must be compiled with parallel functionality");


    static const std::array<hid_t, 12> c_types =
            {0, H5T_NATIVE_CHAR,
                    H5T_NATIVE_UINT8, H5T_NATIVE_INT8,
                    H5T_NATIVE_UINT16, H5T_NATIVE_INT16,
                    H5T_NATIVE_UINT32, H5T_NATIVE_INT32,
                    H5T_NATIVE_UINT64, H5T_NATIVE_INT64,
                    H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE};

    std::string type_name(hid_t type);

    template<typename T=void>
    static constexpr uint_t type_ind() { return 0; }

    template<>
    constexpr uint_t type_ind<char>() { return 1; }

    template<>
    constexpr uint_t type_ind<uint8_t>() { return 2; }

    template<>
    constexpr uint_t type_ind<int8_t>() { return 3; }

    template<>
    constexpr uint_t type_ind<uint16_t>() { return 4; }

    template<>
    constexpr uint_t type_ind<int16_t>() { return 5; }

    template<>
    constexpr uint_t type_ind<uint32_t>() { return 6; }

    template<>
    constexpr uint_t type_ind<int32_t>() { return 7; }

    template<>
    constexpr uint_t type_ind<uint64_t>() { return 8; }

    template<>
    constexpr uint_t type_ind<int64_t>() { return 9; }

    template<>
    constexpr uint_t type_ind<float>() { return 10; }

    template<>
    constexpr uint_t type_ind<double>() { return 11; }

    template<typename T>
    const hid_t &type() {
        typedef arith::comp_t<T> comp_t;
        static_assert(type_ind<comp_t>(), "type has no HDF5 equivalent");
        return c_types[type_ind<comp_t>()];
    }

    static bool types_equal(hid_t t1, hid_t t2) {
        return H5Tequal(t1, t2);
    }

    template<typename T>
    bool types_equal(hid_t t2) {
        return types_equal(type<T>(), t2);
    }

    /**
     * @param h5type
     *  type index
     * @return
     *  size in bytes of the type identified by the given HDF5 type index
     */
    hsize_t type_size(hid_t h5type);

    struct StringType {
        const hid_t m_handle;
        const hsize_t m_nchar; // excluding null terminator

    private:
        static hsize_t size_max(const std::vector<std::string>& vec);

        /*
         * use dummy arg so as not to have same prototype as public ctor in case hsize_t coincides with hid_t
         */
        StringType(hsize_t size, int /*dummy*/);
    public:

        StringType(hid_t handle);

        StringType(const std::string& str);

        StringType(const std::vector<std::string>& str_vec);

        ~StringType();

        operator hid_t() const;
    };
}


#endif //M7_HDF5_TYPE_H
