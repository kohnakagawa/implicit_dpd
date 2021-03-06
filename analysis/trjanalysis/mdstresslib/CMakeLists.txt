cmake_minimum_required(VERSION 2.8)

enable_language(C)
enable_language(CXX)
enable_language(Fortran)

set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)

FILE(GLOB MDSTRESS_SOURCES_H   include/*.h)
FILE(GLOB MDSTRESS_SOURCES_CPP src/*.cpp)
SET(MDSTRESS_SOURCES ${MDSTRESS_SOURCES_CPP} ${MDSTRESS_SOURCES_H})

INCLUDE_DIRECTORIES(include/)

# FIND_LIBRARY(LAPACK_LIBRARY NAMES lapack liblapack HINTS /usr/lib REQUIRED)
FIND_LIBRARY(LAPACK_LIBRARY NAMES mkl_intel_lp64 libmkl_intel_lp64 mkl_sequential libmkl_sequential mkl_core libmkl_core PATHS /cm/shared/apps/intel/composer_xe/2015.5.223/mkl/lib/intel64 REQUIRED)

# ADD_LIBRARY(mdstress SHARED ${MDSTRESS_SOURCES})
ADD_LIBRARY(mdstress STATIC ${MDSTRESS_SOURCES})

TARGET_LINK_LIBRARIES(mdstress ${LAPACK_LIBRARY} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

# INSTALL(TARGETS mdstress RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib/static)

# INSTALL (FILES ${MDSTRESS_SOURCES_H} DESTINATION include/mdstress)

# ADD_SUBDIRECTORY(boostpython)
# ADD_SUBDIRECTORY(tensortools)
