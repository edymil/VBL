cmake_minimum_required(VERSION 3.1)
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
   #this one is needed for the project itself
   #SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/milotti_mts/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# Include all the transitive dependencies determined by Spack.
# If we're not running with Spack, this does nothing...
include_directories($ENV{SPACK_TRANSITIVE_INCLUDE_PATH})

###############  end spack stuff  ##############

# c++ 11 for newer cmake
set(CMAKE_CXX_STANDARD 11)
#always use new c++11 the old way
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

### it is NOT elegant, but we do not know how to handle gmp on mac, 
### we set the compiler flags manually there
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
  if(${APPLE})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG -lgmpxx -lmpfr -lgmp ")
  endif()

endif()

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
find_package(Boost REQUIRED COMPONENTS "thread")
set(myLibs "${myLibs};${Boost_THREAD_LIBRARY}")
message("Boost_THREAD_LIBRARY library: ${Boost_THREAD_LIBRARY}")

### CGAL
find_package(CGAL REQUIRED)
if(CGAL_FOUND)
  message("Found CGAL_LIBRARY: ${CGAL_LIBRARY})")
  set(myLibs "${myLibs};${CGAL_LIBRARY}" )
endif()

#execute_process(COMMAND bash -c "/Users/spack/spack/bin/spack location -i cgal" COMMAND tr -d '\n'  OUTPUT_VARIABLE CGAL_ROOT)
#include("${CGAL_ROOT}/lib/CGAL/FindGMP")


########################################## FindGMP
# This module defines:
#  GMP_FOUND             - system has GMP lib
#  GMP_INCLUDE_DIR       - the GMP include directory
#  GMP_LIBRARIES_DIR     - directory where the GMP libraries are located
#  GMP_LIBRARIES         - Link these to use GMP
#  GMP_IN_CGAL_AUXILIARY - TRUE if the GMP found is the one distributed with CGAL in the auxiliary folder

# TODO: support MacOSX

include(FindPackageHandleStandardArgs)
include(CGAL_GeneratorSpecificSettings)

if(GMP_INCLUDE_DIR)
  set(GMP_in_cache TRUE)
else()
  set(GMP_in_cache FALSE)
endif()
if(NOT GMP_LIBRARIES)
  set(GMP_in_cache FALSE)
endif()

# Is it already configured?
if (GMP_in_cache)

  set(GMP_FOUND TRUE)

else()

  find_path(GMP_INCLUDE_DIR
            NAMES gmp.h
            HINTS ENV GMP_INC_DIR
                  ENV GMP_DIR
                  ${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/include
            PATH_SUFFIXES include
  	        DOC "The directory containing the GMP header files"
           )

  if ( GMP_INCLUDE_DIR STREQUAL "${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/include" )
    cache_set( GMP_IN_CGAL_AUXILIARY TRUE )
  endif()

  find_library(GMP_LIBRARIES NAMES gmp libgmp-10
    HINTS ENV GMP_LIB_DIR
          ENV GMP_DIR
          ${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/lib
    PATH_SUFFIXES lib
    DOC "Path to the GMP library"
    )

  if ( GMP_LIBRARIES )
    get_filename_component(GMP_LIBRARIES_DIR ${GMP_LIBRARIES} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for GMP if couldn't be found
  if ( NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARIES_DIR )
    include( GMPConfig OPTIONAL )
  endif()

  find_package_handle_standard_args(GMP "DEFAULT_MSG" GMP_LIBRARIES GMP_INCLUDE_DIR)

endif()
########################################## FindGMP



### MPFR
##set(myLibs "${myLibs};mpfr")
## find_package(MPFR REQUIRED)
#find_library(MPFR_LIBRARY gmp)
##if(MPFR_FOUND)
#  message("Found MPFR_LIBRARY: ${MPFR_LIBRARY})")
#  set(myLibs "${myLibs};${MPFR_LIBRARY}" )
##else()
#  message("MPFR library: ${GMP_LIBRARY}")
##  message(FATAL_ERROR "MPFR stuff is not found, but labeled as found")
##endif()


### cmake debug output CMAKE_EXE_LINKER_FLAGS
message("myLibs: ${myLibs}")
message("compiler flags: ${CMAKE_CXX_FLAGS}")
message("CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")


add_subdirectory(src)
