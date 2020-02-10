#include <iostream>
#include <iomanip>
#include "data/DataTable.h"
#include "data/BitfieldNew.h"
#include "parallel/MPIWrapper.h"
#include "datasystems/DataSystem.h"
#include "io/InputOptions.h"
#include "CLI/CLI.hpp"

int main(int argc, char **argv) {

    MPIWrapper::initialize(&argc, &argv);
    MPIWrapper mpi;

    /*
     * Setup and read-in runtime options from the command line
     */
    CLI::App cli_app{InputOptions::description};
    InputOptions input(cli_app);
    if(!mpi.i_am_root()){
        /*
         * only allow standard output and error from the root MPI rank
         */
        std::cout.setstate(std::ios_base::failbit);
        std::cerr.setstate(std::ios_base::failbit);
    }
    try {
        cli_app.parse((argc), (argv));                                                                                   \

    } catch (const CLI::ParseError &e) {
        MPIWrapper::finalize();
        cli_app.exit(e);
        return 0;
    }

    std::array<size_t, nnumeric> numeric_lengths;
    //std::vector<size_t> bitfield_lengths{12, 18};
    numeric_lengths[type_number<int>] = 3;
    //numeric_lengths[index<std::complex<float>>] = 2;
    numeric_lengths[type_number<float>] = 2;

    //DataTable dataTable(numeric_lengths, bitfield_lengths, 5, 2);
    DataTable store(numeric_lengths, defs::inds{}, 5, 1);
    DataTable send(numeric_lengths, defs::inds{}, 5, mpi.nrank());
    DataTable recv(numeric_lengths, defs::inds{}, 5, 1);

    DataSystem dataSystem(store, send, recv);

    size_t irow;

    if (mpi.i_am(0)) {
        irow = send.claim_rows(0);
        *send.view<int>(0, irow, 0) = 9000;
        *send.view<int>(0, irow, 1) = 9001;
        *send.view<int>(0, irow, 2) = 9002;

        irow = send.claim_rows(1);
        *send.view<int>(1, irow, 0) = 9010;
        *send.view<int>(1, irow, 1) = 9011;
        *send.view<int>(1, irow, 2) = 9012;

        irow = send.claim_rows(1);
        *send.view<int>(1, irow, 0) = 9013;
        *send.view<int>(1, irow, 1) = 9014;
        *send.view<int>(1, irow, 2) = 9015;
    }

    if (mpi.i_am(1)) {
        irow = send.claim_rows(0);
        *send.view<int>(0, irow, 0) = 9100;
        *send.view<int>(0, irow, 1) = 9101;
        *send.view<int>(0, irow, 2) = 9102;
    }

    for (auto i{0ul}; i<mpi.nrank(); ++i){
        mpi.barrier();
        if (mpi.i_am(i)) send.print();
    }

    dataSystem.communicate();

    for (auto i{0ul}; i<mpi.nrank(); ++i){
        mpi.barrier();
        if (mpi.i_am(i)) recv.print();
    }

    MPIWrapper::finalize();

    return 0;

}