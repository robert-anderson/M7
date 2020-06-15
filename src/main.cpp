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
    CLI::App cli_app{InputOptions::description};
    InputOptions input(cli_app);


    std::ofstream ofstdout, ofstderr;
    if (!mpi::i_am_root()) {
        /*
         * only allow standard output and error from the root MPI rank
         */
#ifdef DNDEBUG
        std::cout.setstate(std::ios_base::failbit);
        std::cerr.setstate(std::ios_base::failbit);
#else
        ofstdout = std::ofstream("rank_"+std::to_string(mpi::irank())+".out");
        std::cout.rdbuf(ofstdout.rdbuf());
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

    FciqmcCalculation calc(input);
    calc.execute();

    mpi::finalize();

    return 0;

}
