#include <iostream>
#include <iomanip>
#include <M7_lib/parallel/MPIWrapper.h>
#include <M7_lib/dynamics/FciqmcCalculation.h>
#include <M7_lib/field/Mbf.h>
#include <M7_lib/conf/Conf.h>

int main(int argc, char **argv) {

    mpi::initialize(&argc, &argv);

    logging::title();

    if (argc == 1) {
        // input file not provided, print out help string
        for (auto line: logging::make_defs_table()) std::cout << line << std::endl;
        conf::Document().print_help();
        mpi::finalize();
        return 0;
    }

    logging::initialize();
    /*
     * log compile-time defintions
     */
    logging::info_lines(logging::make_defs_table());

    logging::info("Number of MPI ranks in world communicator: {}", mpi::g_nrank);
    logging::info("Number of MPI ranks per node: {}", mpi::g_nrank_on_node);

    conf::Document opts(argv[1]);
    opts.validate();
    opts.log();

    std::streambuf *original_stdout_buffer = nullptr;
    std::streambuf *original_stderr_buffer = nullptr;
    if (!mpi::i_am_root()) {
        /*
         * only allow standard output and error from the root MPI rank
         */
        std::cout.setstate(std::ios_base::failbit);
        std::cerr.setstate(std::ios_base::failbit);
    }

    {
        mpi::barrier();
        fciqmc::run(opts);
    }

    if (!mpi::i_am_root()) {
        /*
         * reinstate original buffers
         */
        std::cout.rdbuf(original_stdout_buffer);
        std::cerr.rdbuf(original_stderr_buffer);
    }
    mpi::finalize();

    return 0;

}
