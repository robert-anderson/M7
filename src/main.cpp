#include <iostream>
#include <iomanip>
#include <src/core/parallel/MPIWrapper.h>
#include <src/core/dynamics/FciqmcCalculation.h>
#include <src/core/io/InputOptions.h>
#include "CLI/CLI.hpp"

int main(int argc, char **argv) {

    mpi::initialize(&argc, &argv);
    /*
     * Setup and read-in runtime options from the command line
     */
    CLI::App cli_app{InputOptions::program_description};
    InputOptions input(cli_app);


    std::streambuf *original_stdout_buffer = nullptr;
    std::streambuf *original_stderr_buffer = nullptr;
    std::ofstream ofstdout, ofstderr;
    if (!mpi::i_am_root()) {
        /*
         * only allow standard output and error from the root MPI rank
         */
#ifdef DNDEBUG
        std::cout.setstate(std::ios_base::failbit);
        std::cerr.setstate(std::ios_base::failbit);
#else
        original_stdout_buffer = std::cout.rdbuf();
        ofstdout = std::ofstream("rank_"+std::to_string(mpi::irank())+".out");
        std::cout.rdbuf(ofstdout.rdbuf());
        original_stderr_buffer = std::cerr.rdbuf();
        ofstderr = std::ofstream("rank_"+std::to_string(mpi::irank())+".err");
        std::cerr.rdbuf(ofstderr.rdbuf());
#endif
    }

    try {
        cli_app.parse((argc), (argv));
    } catch (const CLI::ParseError &e) {
        mpi::finalize();
        cli_app.exit(e);
        return 0;
    }

    {
        FciqmcCalculation calc(input);
        calc.execute();
    }

    if (!mpi::i_am_root()) {
        /*
         * reinstate original buffers
         */
#ifndef DNDEBUG
        std::cout.rdbuf(original_stdout_buffer);
        std::cerr.rdbuf(original_stderr_buffer);
#endif
    }
    mpi::finalize();

    return 0;

}
