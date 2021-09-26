#include <iostream>
#include <iomanip>
#include <src/core/parallel/MPIWrapper.h>
#include <src/core/dynamics/FciqmcCalculation.h>
#include <src/core/io/YamlWrapper.h>
#include <src/core/config/FciqmcConfig.h>

int main(int argc, char **argv) {

    mpi::initialize(&argc, &argv);
    log::initialize();

    log::info("Initializing M7 instance");
    /*
     * log compile-time defintions
     */
    log::info("Many-body basis definition: {}", consts::mbf_type_name<defs::mbf_ind>());
    log::info("Walker arithmetic type: {}", defs::enable_complex ? "complex" : "real");

    if (argc == 1) {
        // input file not provided, print out help string
        std::cout << fciqmc_config::Document(nullptr).help_string() << std::endl;
        mpi::finalize();
        return 0;
    }

    auto yf = yaml::File(std::string(argv[1]));
    fciqmc_config::Document opts(&yf);
    opts.verify();
    opts.log_value();

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
