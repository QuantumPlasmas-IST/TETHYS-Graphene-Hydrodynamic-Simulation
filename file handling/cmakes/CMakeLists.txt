cmake_minimum_required (VERSION 3.10.1)
project( myFirstHdf5 C CXX )

set (LIB_TYPE STATIC) # or SHARED
string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)

#find_package (HDF5 NAMES hdf5 COMPONENTS C CXX ${SEARCH_TYPE})
#find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)

link_directories( ${HDF5_LIBRARY_DIRS} )

include_directories (${HDF5_INCLUDE_DIR})
set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

#set (example hdfcompile)

add_executable (myFirstHdf5 myFirstHdf5.cpp)

target_link_libraries (myFirstHdf5 ${HDF5_CXX_LIBRARIES})
