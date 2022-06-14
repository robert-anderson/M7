target_sources(unittest
    PRIVATE
        BosonCouplings.cpp
        TranscorrelatedHamiltonian.cpp
        TranscorrelatedBosonHamiltonian.cpp
)

add_subdirectory(bos)
add_subdirectory(frm)
add_subdirectory(frmbos)
