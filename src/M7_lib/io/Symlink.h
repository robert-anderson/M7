//
// Created by rja on 27/05/22.
//

#ifndef M7_SYMLINK_H
#define M7_SYMLINK_H

#include "M7_lib/defs.h"

/**
 * simple wrapper to the UNIX standard lib to create a soft link to a file located at "src", to a given path "dst"
 * dtor calls unlink so the instance cleans up after itself when it goes out of scope.
 * this class is useful in unit tests of interfaces to external programs which require inputs to be available in the
 * current working directory, but working directory-agnostic code is *strongly* preferred elsewhere
 */
struct Symlink {
    const str_t m_fname;
    Symlink(str_t src, str_t dst);
    ~Symlink();
};

/**
 * convenient specialization to the commonly-encountered case in testing where the src is specified relative to the
 * "assets" directory
 */
struct AssetSymlink : Symlink {
    AssetSymlink(str_t src, str_t dst);
};

#endif //M7_SYMLINK_H
