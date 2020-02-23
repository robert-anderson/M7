//
// Created by Robert John Anderson on 2020-02-23.
//

#include "StatsFile.h"

void StatsFile::write_header() {
    if (mpi::i_am_root()) {
        size_t ncolumn = 0ul;
        for (auto &column : m_columns) ncolumn += column.nsubcolumn();
        m_file <<
               "################################\n"
               "#    M7 FCIQMC Stats Output    #\n"
               "################################\n"
               "#\n"
               "# Number of statistics output: " << m_columns.size() << "\n"
                                                                        "# Number of columns: " << ncolumn
               << "\n#\n";
        size_t icol = 1ul;
        for (auto &column : m_columns) {
            if (numtypes::is_complex(column.itype)) {
                m_file << "#  " << icol++ << ".  " << column.description << " (real)\n";
                m_file << "#  " << icol++ << ".  " << column.description << " (imag)\n";
            } else {
                m_file << "#  " << icol++ << ".  " << column.description << " \n";
            }
        }
        m_file << std::endl;
    }
}

void StatsFile::flush() {
    if (mpi::i_am_root()) {
        for (auto &column : m_columns) {
            m_file << column.output << "  ";
            column.output = "none";
        }
        m_file << std::endl;
    }
}
