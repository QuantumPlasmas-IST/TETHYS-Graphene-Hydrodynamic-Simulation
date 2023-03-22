cmake_minimum_required(VERSION 3.16)
project(TETHYS_Hydrodynamic_Simulation)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")

enable_testing()

include_directories(src/)
include_directories(src/includes)

find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})

#find_package(GSL REQUIRED)
#link_libraries(GSL::gsl GSL::gslcblas)

#set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -O2")
set(CMAKE_CXX_FLAGS "-O2")
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS}")
endif()



add_executable(TETHYS_1D src/TETHYS_1D_Main_v262.cpp)
add_executable(TETHYS_2D src/TETHYS_2D_Main_v262.cpp)
add_executable(TETHYS_2D_CNP src/TETHYS_2D_Dirac_main.cpp)
add_executable(TETHYS_ELEC_2D src/TETHYS_2D_ElectronicAnalysis.cpp)
#add_executable(TETHYS_RADIATION_2D src/TETHYS_2D_Radiation.cpp)

file(GLOB SRC_LIB_FILES src/*Lib.cpp)



add_library(tethys_library STATIC ${SRC_LIB_FILES})


target_link_libraries(TETHYS_1D PUBLIC tethys_library)

target_link_libraries(TETHYS_2D PUBLIC tethys_library)
target_link_libraries(TETHYS_2D_CNP PUBLIC tethys_library)
target_link_libraries(TETHYS_ELEC_2D PUBLIC tethys_library)

#target_link_libraries(TETHYS_RADIATION_2D PUBLIC tethys_library)


target_include_directories(TETHYS_1D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_1D PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TETHYS_2D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_2D PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TETHYS_2D_CNP PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_2D_CNP PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TETHYS_ELEC_2D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_ELEC_2D PRIVATE ${HDF5_CXX_LIBRARIES})


add_test(test1D TETHYS_1D 23 11 0.8 0.04 0 0 0 1 1)
add_test(test2D TETHYS_2D 2DTESTparameters.ini)
set(H5_TEST_FILE "hdf5_2D_S=31.00vF=10.50vis=0.020odd=0.010l=0.400wc=0.70therm=0.00.h5")
add_test(testElec2D TETHYS_ELEC_2D ${H5_TEST_FILE})
add_test(test2DCNP TETHYS_2D_CNP 2DTESTCNPparameters.ini)





# Doxygen

# look for Doxygen package
find_package(Doxygen)

if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.out)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # Note: do not put "ALL" - this builds docs together with application EVERY TIME!
    add_custom_target( docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating documentation with Doxygen"
            VERBATIM )
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)