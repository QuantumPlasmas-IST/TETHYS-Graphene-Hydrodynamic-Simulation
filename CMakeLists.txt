cmake_minimum_required(VERSION 3.16)
project(Hydrodynamic_Simulation)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

find_package(HDF5 REQUIRED COMPONENTS CXX)

find_package(GSL REQUIRED)
link_libraries(GSL::gsl GSL::gslcblas)


set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -O2")
#set(GSL_FLAGS "-lgsl -lgslcblas -lm")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=leak -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined -g ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory -g ")

add_executable(TETHYS_1D
        BoundaryLib.cpp
        BoundaryLib.h
        ElectricLib.cpp
        ElectricLib.h
        Tethys1DLib.cpp
        Tethys1DLib.h
        Tethys2DLib.cpp
        Tethys2DLib.h
        TETHYS_1D_Main_v134.cpp
        TethysLib.cpp
        TethysLib.h)

add_executable(TETHYS_2D
        BoundaryLib.cpp
        BoundaryLib.h
        ElectricLib.cpp
        ElectricLib.h
        Tethys1DLib.cpp
        Tethys1DLib.h
        Tethys2DLib.cpp
        Tethys2DLib.h
        TETHYS_2D_Main_v210.cpp
        TethysLib.cpp
        TethysLib.h)





target_include_directories(TETHYS_1D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_1D PRIVATE ${HDF5_CXX_LIBRARIES})

target_include_directories(TETHYS_2D PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(TETHYS_2D PRIVATE ${HDF5_CXX_LIBRARIES})
