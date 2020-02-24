#include <iostream>
#include <iomanip>
#include <src/data/BitfieldNew.h>
#include <src/parallel/MPIWrapper.h>
#include <src/dynamics/FciqmcCalculation.h>
#include <src/io/InputOptions.h>
#include "CLI/CLI.hpp"


#include <src/enumerators/VectorCombinationEnumerator.h>

int main(int argc, char **argv) {

    mpi::initialize(&argc, &argv);
    /*
     * Setup and read-in runtime options from the command line
     */
    CLI::App cli_app{InputOptions::description};
    InputOptions input(cli_app);
    if (!mpi::i_am_root()) {
        /*
         * only allow standard output and error from the root MPI rank
         */
        std::cout.setstate(std::ios_base::failbit);
        std::cerr.setstate(std::ios_base::failbit);
    }
    try {
        cli_app.parse((argc), (argv));                                                                                   \

    } catch (const CLI::ParseError &e) {
        mpi::finalize();
        cli_app.exit(e);
        return 0;
    }

    FciqmcCalculation calc(input);
    calc.run();

    mpi::finalize();

    return 0;

}