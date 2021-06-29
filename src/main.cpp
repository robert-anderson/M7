#include <iostream>
#include <iomanip>
#include <src/core/parallel/MPIWrapper.h>
#include <src/core/dynamics/FciqmcCalculation.h>
#include <src/core/io/InputOptions.h>
#include <src/core/io/YamlWrapper.h>
#include <src/core/config/FciqmcConfig.h>
#include "CLI/CLI.hpp"

int main(int argc, char **argv) {

    mpi::initialize(&argc, &argv);
    log::initialize();
    /*
     * Setup and read-in runtime options from the command line
     */
    CLI::App cli_app{InputOptions::program_description};
    InputOptions input(cli_app);

    if (argc == 1){
        // input file not provided, print out help string
        std::cout << fciqmc_config::Document(nullptr).help_string() << std::endl;
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

    try {
        cli_app.parse((argc), (argv));
    } catch (const CLI::ParseError &e) {
        mpi::finalize();
        cli_app.exit(e);
        return 0;
    }

    {
        input.init();
        FciqmcCalculation calc(input);
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
