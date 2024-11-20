cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set( PROJECTNAME "gebt" )

enable_language( Fortran )

project( ${PROJECTNAME} Fortran )


# set( CMAKE_CXX_STANDARD 11 )

set( VERSION_MAJOR 4 )
set( VERSION_MINOR 8 )
set( VERSION_FULL "${VERSION_MAJOR}.${VERSION_MINOR}" )



# Source files
# ------------
set( PROJ_SOURCES 
  src/ddep.f
  src/ma28.f
  src/mc19.f
  # src/blas.f
  # src/lapack.f
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

set(EXTRA_LIB)


# BLAS and LAPACK
# ---------------
# find_library( BLAS NAMES blas REQUIRED )
# message( ${BLAS} )
# find_library( LAPACK NAMES lapack REQUIRED )
# message( ${LAPACK} )

find_package( BLAS )
if (BLAS_FOUND)
  message( STATUS "BLAS found: ${BLAS_LIBRARIES}" )
  list(APPEND EXTRA_LIB ${BLAS_LIBRARIES})
else()
  message( STATUS "BLAS not found, using src/blas.f" )
  list(APPEND PROJ_SOURCES src/blas.f)
endif()

find_package( LAPACK )
if (LAPACK_FOUND)
  message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}" )
  list(APPEND EXTRA_LIB ${LAPACK_LIBRARIES})
else()
  message( STATUS "LAPACK not found, using src/lapack.f" )
  list(APPEND PROJ_SOURCES src/lapack.f)
endif()


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


add_executable( ${PROJECTNAME} ${PROJ_SOURCES} )
# add_executable( ${PROJECTNAME} ${MA28_SOURCES} ${arpackutil_STAT_SRCS} ${arpacksrc_STAT_SRCS} ${PROJ_SOURCES} )

if (EXTRA_LIB)
  target_link_libraries( ${PROJECTNAME} ${EXTRA_LIB} )
endif()

# Install
# -------
install(TARGETS ${PROJECTNAME} RUNTIME DESTINATION bin)
