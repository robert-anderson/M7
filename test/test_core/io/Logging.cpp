//
// Created by anderson on 18/07/2022.
//

#include <M7_lib/io/Logging.h>
#include <M7_lib/io/FileReader.h>
#include "gtest/gtest.h"

#ifdef ENABLE_LOCAL_LOGGING
TEST(Logging, Local) {
    logging::info_("i am rank {}", mpi::irank());
    logging::error_("i am rank {}", mpi::irank());
    logging::warn_("i am rank {}", mpi::irank());
    logging::critical_("i am rank {}", mpi::irank());
    logging::flush_all();
    /*
    [TIMESTAMP] [info] i am rank 0
    [TIMESTAMP] [error] i am rank 0
    [TIMESTAMP] [warning] i am rank 0
    [TIMESTAMP] [critical] i am rank 0
     */

    FileReader reader(mpi::nrank()==1ul ? str_t("M7.log") : logging::format("M7.log.{}", mpi::irank()));
    str_t line;
    reader.next(line);
    ASSERT_NE(line.find("[info]"), ~0ul);
    reader.next(line);
    ASSERT_NE(line.find("[error]"), ~0ul);
    reader.next(line);
    ASSERT_NE(line.find("[warning]"), ~0ul);
    reader.next(line);
    ASSERT_NE(line.find("[critical]"), ~0ul);
}
#endif