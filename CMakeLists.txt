cmake_minimum_required(VERSION 2.8)
project(particleBasedTumor)
set(CMAKE_CXX_STANDARD 11) # C++11...
### compiler
message("CMAKE_CXX_COMPILER:        ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_COMPILER_ID:     ${CMAKE_CXX_COMPILER_ID}")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  message("using GNU build system")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math")
endif()
### Boost
find_package(Boost REQUIRED COMPONENTS "thread")
set(myLibs "${myLibs};${Boost_THREAD_LIBRARY}")

### OPENMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  message(STATUS "OpenMP_CXX_FLAGS      = ${OpenMP_CXX_FLAGS}")
  message(STATUS "OpenMP_C_FLAGS      = ${OpenMP_C_FLAGS}")
  set (CMAKE_C_FLAGS  "${OpenMP_C_FLAGS}   -DOPENMP_ENABLED=1 ${CMAKE_C_FLAGS}  ")
  set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -DOPENMP_ENABLED=1 ${CMAKE_CXX_FLAGS}")
endif()
### CGAL
find_package(CGAL REQUIRED)
if(CGAL_FOUND)
  message("Found CGAL_LIBRARY: ${CGAL_LIBRARY})")
  set(myLibs "${myLibs};${CGAL_LIBRARY}" )
endif()

### GMP
set(myLibs "${myLibs};gmp")

message("set comiler flags: ${CMAKE_CXX_FLAGS}")

add_subdirectory(src)
