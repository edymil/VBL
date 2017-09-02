cmake_minimum_required(VERSION 2.8)
project(VBL)

###############  set up spack stuff  ##############

# enable @rpath in the install name for any shared library being built
# note: it is planned that a future version of CMake will enable this by default
set(CMAKE_MACOSX_RPATH 1)

# Always use full RPATH
# http://www.cmake.org/Wiki/CMake_RPATH_handling
# http://www.kitware.com/blog/home/post/510

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# Include all the transitive dependencies determined by Spack.
# If we're not running with Spack, this does nothing...
include_directories($ENV{SPACK_TRANSITIVE_INCLUDE_PATH})

###############  end spack stuff  ##############

# c++ 11 for newer cmake
set(CMAKE_CXX_STANDARD 11)
#always use new c++11 the old way
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

### check if apple
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(APPLE TRUE)
    message("Your are using an apple PC.")
endif()

### compiler
message("CMAKE_CXX_COMPILER:        ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_COMPILER_ID:     ${CMAKE_CXX_COMPILER_ID}")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  message("using GNU build system")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG ")

set(myIncs "")
set(myLibs "")

### OPENMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  message(STATUS "OpenMP_CXX_FLAGS        = ${OpenMP_CXX_FLAGS}")
  message(STATUS "OpenMP_C_FLAGS          = ${OpenMP_C_FLAGS}")
  message(STATUS "OpenMP_EXE_LINKER_FLAGS = ${OpenMP_EXE_LINKER_FLAGS}")
  set (CMAKE_C_FLAGS          "${OpenMP_C_FLAGS} ${CMAKE_C_FLAGS}")
  set (CMAKE_CXX_FLAGS        "${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  #set (CMAKE_C_FLAGS  "${OpenMP_C_FLAGS}   -DOPENMP_ENABLED=1 ${CMAKE_C_FLAGS}  ")
  #set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -DOPENMP_ENABLED=1 ${CMAKE_CXX_FLAGS}")
  #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  message(FATAL_ERROR "OpenMP stuff is not found, but labeled as found")
endif()

### Boost
execute_process(COMMAND bash -c "/Users/spack/spack/bin/spack location -i boost" COMMAND tr -d '\n'  OUTPUT_VARIABLE BOOST_ROOT)
message("Boost directory ${BOOST_ROOT}")
set(LIBS "${BOOST_ROOT}/lib")
set(INCS "${BOOST_ROOT}/include")

### GMP
execute_process(COMMAND bash -c "/Users/spack/spack/bin/spack location -i gmp" COMMAND tr -d '\n'  OUTPUT_VARIABLE GMP_ROOT)
message("GMP directory ${GMP_ROOT}")
set(LIBS "${LIBS};${GMP_ROOT}/lib")
set(INCS "${INCS};${GMP_ROOT}/include")

### MPFR
execute_process(COMMAND bash -c "/Users/spack/spack/bin/spack location -i mpfr" COMMAND tr -d '\n'  OUTPUT_VARIABLE MPFR_ROOT)
message("MPFR directory ${MPFR_ROOT}")
set(LIBS "${LIBS};${MPFR_ROOT}/lib")
set(INCS "${INCS};${MPFR_ROOT}/include")

### CGAL
execute_process(COMMAND bash -c "/Users/spack/spack/bin/spack location -i cgal" COMMAND tr -d '\n'  OUTPUT_VARIABLE CGAL_ROOT)
message("CGAL directory ${CGAL_ROOT}")
set(LIBS "${LIBS};${CGAL_ROOT}/lib")
set(INCS "${INCS};${CGAL_ROOT}/include")

### cmake debug output CMAKE_EXE_LINKER_FLAGS
message("set content of mylib: ${myLibs}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${INCS} -L${LIBS} -lboost-thread-mt -lgmp -lgmpxx -lmpfr")
message("CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message("CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")


add_subdirectory(src)

