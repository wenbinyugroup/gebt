cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set( PROJECTNAME "gebt" )

project( ${PROJECTNAME} Fortran )


# set( CMAKE_CXX_STANDARD 11 )

set( VERSION_MAJOR 4 )
set( VERSION_MINOR 8 )
set( VERSION_FULL "${VERSION_MAJOR}.${VERSION_MINOR}" )

# include(CheckCXXCompilerFlag)
# CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
# CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
# if(COMPILER_SUPPORTS_CXX11)
# 	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# elseif(COMPILER_SUPPORTS_CXX0X)
# 	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
# else()
#   message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
# endif()


# BLAS and LAPACK
# ---------------
# find_library( BLAS NAMES blas REQUIRED )
# message( ${BLAS} )
# find_library( LAPACK NAMES lapack REQUIRED )
# message( ${LAPACK} )

find_package( LAPACK REQUIRED )
message( STATUS ${LAPACK_INCLUDE_DIRS} )
message( STATUS ${LAPACK_LIBRARY_DIRS} )
message( STATUS ${LAPACK_LIBRARIES} )

find_package( BLAS REQUIRED )
message( STATUS ${BLAS_INCLUDE_DIRS} )
message( STATUS ${BLAS_LIBRARY_DIRS} )
message( STATUS ${BLAS_LIBRARIES} )


# MA28 (linear solver)
# --------------------
# set( MA28_SOURCES 
#   src/ddep.f
#   src/ma28.f
#   src/mc19.f
# )


# ARPACK (eigen solver)
# ---------------------
# set(arpack_SOURCE_DIR arpack-ng)

# include_directories( ${CMAKE_SOURCE_DIR}/include )
# file(GLOB arpacksrc_STAT_SRCS ${arpack_SOURCE_DIR}/dbgini.f ${arpack_SOURCE_DIR}/staini.f ${arpack_SOURCE_DIR}/SRC/*.f)

# set(arpackutil_STAT_SRCS
#     ${arpack_SOURCE_DIR}/UTIL/icnteq.f
#     ${arpack_SOURCE_DIR}/UTIL/icopy.f
#     ${arpack_SOURCE_DIR}/UTIL/iset.f
#     ${arpack_SOURCE_DIR}/UTIL/iswap.f
#     ${arpack_SOURCE_DIR}/UTIL/ivout.f
#     ${arpack_SOURCE_DIR}/UTIL/second_NONE.f
#     ${arpack_SOURCE_DIR}/UTIL/svout.f
#     ${arpack_SOURCE_DIR}/UTIL/smout.f
#     ${arpack_SOURCE_DIR}/UTIL/dvout.f
#     ${arpack_SOURCE_DIR}/UTIL/dmout.f
#     ${arpack_SOURCE_DIR}/UTIL/cvout.f
#     ${arpack_SOURCE_DIR}/UTIL/cmout.f
#     ${arpack_SOURCE_DIR}/UTIL/zvout.f
#     ${arpack_SOURCE_DIR}/UTIL/zmout.f )


# GEBT
# ----
set( PROJ_SOURCES 
  src/ddep.f
  src/ma28.f
  src/mc19.f
  src/arpack.f
  src/CPUtime.f90
  src/GlobalDataFun.f90
  src/TimeFunction.f90
  src/PrescribedCondition.f90
  src/InternalData.f90
  src/Preprocess.f90
  src/Element.f90
  src/Member.f90
  src/System.f90
  src/Solve.f90
  src/EigenSolve.f90
  src/Analysis.f90
  src/IO.f90
  src/main.f90
)

add_executable( ${PROJECTNAME} ${PROJ_SOURCES} )
# add_executable( ${PROJECTNAME} ${MA28_SOURCES} ${arpackutil_STAT_SRCS} ${arpacksrc_STAT_SRCS} ${PROJ_SOURCES} )

target_link_libraries( ${PROJECTNAME} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )

# Install
# -------
install(TARGETS ${PROJECTNAME} RUNTIME DESTINATION bin)
