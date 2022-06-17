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


#add_executable(TEST_laplacian src/laplacian_main.cpp)


add_executable(TETHYS_1D src/TETHYS_1D_Main_v134.cpp)
add_executable(TETHYS_2D src/TETHYS_2D_Main_v261.cpp)
add_executable(TETHYS_ELEC_2D src/TETHYS_2D_ElectronicAnalysis.cpp)
add_executable(TETHYS_RADIATION_2D src/TETHYS_2D_Radiation.cpp)

#file(GLOB SRC_LIB_FILES src/*Lib.cpp)
file(GLOB SRC_LIB_FILES src/*Lib.cpp src/Library_Radiation/*.cpp)

add_library(tethys_library STATIC ${SRC_LIB_FILES})

#target_link_libraries(TEST_laplacian PUBLIC tethys_library)

target_link_libraries(TETHYS_1D PUBLIC tethys_library)

target_link_libraries(TETHYS_2D PUBLIC tethys_library)
target_link_libraries(TETHYS_ELEC_2D PUBLIC tethys_library)

target_link_libraries(TETHYS_RADIATION_2D PUBLIC tethys_library)

#target_include_directories(TEST_laplacian PRIVATE ${HDF5_INCLUDE_DIRS})
#target_link_libraries(TEST_laplacian PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TETHYS_1D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_1D PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TETHYS_2D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_2D PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TETHYS_ELEC_2D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_ELEC_2D PRIVATE ${HDF5_CXX_LIBRARIES})

add_test(test2D TETHYS_2D 2DTESTparameters.ini)
set(H5_TEST_FILE "hdf5_2D_S=21.00vF=10.50vis=0.010odd=0.002l=0.001wc=0.00therm=0.45.h5")
add_test(testElec2D TETHYS_ELEC_2D ${H5_TEST_FILE})