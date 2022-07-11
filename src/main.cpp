#include <iostream>
#include <iomanip>
#include <M7_lib/parallel/MPIWrapper.h>
#include <M7_lib/dynamics/FciqmcCalculation.h>
#include <M7_lib/field/Mbf.h>
#include <M7_lib/conf/Conf.h>

int main(int argc, char **argv) {

    mpi::initialize(&argc, &argv);
    log::initialize();

    log::info("Initializing M7 instance");
    /*
     * log compile-time defintions
     */
    log::info("M7 {} build", c_enable_debug ? "DEBUG" : "RELEASE");
    log::info("Many-body basis definition: {}", mbf::name<c_mbf_type_ind>());
    log::info("Walker arithmetic type: {}", c_enable_complex_wf ? "complex" : "real");

    if (argc == 1) {
        // input file not provided, print out help string
        conf::Document().print_help();
        mpi::finalize();
        return 0;
    }

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
        FciqmcCalculation calc(opts);
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
