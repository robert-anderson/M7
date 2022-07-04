//
// Created by rja on 27/05/22.
//

#include "Symlink.h"
#include "unistd.h"
#include "M7_lib/parallel/MPIAssert.h"

Symlink::Symlink(str_t src, str_t dst) : m_fname(std::move(dst)) {
    REQUIRE_EQ(symlink(src.c_str(), m_fname.c_str()), 0, "Error creating symbolic link");
}

Symlink::~Symlink() {
    REQUIRE_EQ(unlink(m_fname.c_str()), 0, "Error deleting symbolic link");
}

AssetSymlink::AssetSymlink(str_t src, str_t dst) : Symlink(PROJECT_ROOT"/assets/" + src, dst){}
