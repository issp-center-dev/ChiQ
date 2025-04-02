cmake_minimum_required(VERSION 3.5...3.16)
project(ChiQ)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(Eigen3 3.1 REQUIRED)

find_package(Threads)

set(CHIQ_LIBRARY_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
include_directories(${CHIQ_LIBRARY_INCLUDE_DIRS})
include_directories(./src)

# Build and enable tests
# testing setup
# enable_testing() must be called in the top-level CMakeLists.txt before any add_subdirectory() is called.
option(Testing "Enable testing" OFF)

if (Testing)
  set(LINK_ALL ${CMAKE_THREAD_LIBS_INIT})
  enable_testing(tests)
  add_subdirectory(tests)
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)

add_subdirectory(extern/pybind11)  # require cmake 3.4
add_subdirectory(src)
add_subdirectory(python)

# Configuration files for ChiQ
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/chiqvars.sh.build.in ${CMAKE_CURRENT_BINARY_DIR}/chiqvars.sh @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/chiqvars.sh.installed.in ${CMAKE_CURRENT_BINARY_DIR}/share/chiqvars.sh @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/share/chiqvars.sh DESTINATION share)
