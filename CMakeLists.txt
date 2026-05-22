cmake_minimum_required(VERSION 3.10)

set( PROJECTNAME "gebt" )

project( ${PROJECTNAME} Fortran )

set( VERSION_MAJOR 4 )
set( VERSION_MINOR 8 )
set( VERSION_FULL "${VERSION_MAJOR}.${VERSION_MINOR}" )

# Source files
# ------------
# MA28/MC19/ddep: Harwell sparse solver — no standard system package.
set( PROJ_SOURCES
  src/ddep.f
  src/ma28.f
  src/mc19.f
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
set(BLAS_LAPACK_LIB)

# Prefer static libraries (.a) over import libraries (.dll.a) so that
# -static linking produces a self-contained executable.
if (WIN32)
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a .dll.a)
endif()

# BLAS and LAPACK
# ---------------
# Try find_package first (works for shared libs and Linux system BLAS).
# On Windows/MinGW with static OpenBLAS, find_package may fail; fall back
# to find_library for libopenblas (which provides both BLAS and LAPACK).
find_package( BLAS )
find_package( LAPACK )
if (BLAS_FOUND AND LAPACK_FOUND)
  message( STATUS "BLAS found: ${BLAS_LIBRARIES}" )
  message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}" )
  set(BLAS_LAPACK_LIB ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
  find_library( OPENBLAS_LIBRARY NAMES openblas )
  if (OPENBLAS_LIBRARY)
    message( STATUS "BLAS/LAPACK: using OpenBLAS ${OPENBLAS_LIBRARY}" )
    set(BLAS_LAPACK_LIB ${OPENBLAS_LIBRARY})
  else()
    message( STATUS "BLAS/LAPACK not found, using bundled src/archive/blas.f + src/archive/lapack.f" )
    list(APPEND PROJ_SOURCES src/archive/blas.f src/archive/lapack.f)
  endif()
endif()

# ARPACK (eigenvalue solver)
# --------------------------
# Prefer system arpack-ng library.  Fall back to bundled src/archive/arpack.f.
find_library( ARPACK_LIBRARY NAMES arpack )
if (ARPACK_LIBRARY)
  message( STATUS "ARPACK: using system library ${ARPACK_LIBRARY}" )
  list(APPEND EXTRA_LIB ${ARPACK_LIBRARY})
else()
  message( STATUS "ARPACK not found, using bundled src/archive/arpack.f" )
  list(APPEND PROJ_SOURCES src/archive/arpack.f)
endif()


add_executable( ${PROJECTNAME} ${PROJ_SOURCES} )

# Compiler-specific flags (gfortran)
# -----------------------------------
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Bundled MA28/MC19/ddep use legacy Fortran constructs; bundled arpack.f
  # (if used) has rank mismatches.  Only -fallow-argument-mismatch is needed
  # when arpack.f is compiled; harmless otherwise.
  target_compile_options(${PROJECTNAME} PRIVATE -fallow-argument-mismatch -std=legacy)
  # Static-link gfortran/libgcc runtimes to avoid a UCRT heap conflict
  # that corrupts TRIM() temporaries (observed with gfortran 16 / MinGW-w64).
  # Do NOT use full -static: it conflicts with OpenMP in OpenBLAS.
  target_link_options(${PROJECTNAME} PRIVATE -static-libgfortran -static-libgcc)
  if (WIN32)
    target_link_options(${PROJECTNAME} PRIVATE -Wl,--stack,67108864)
  endif()
endif()

if (EXTRA_LIB)
  # Link order matters for static libs: arpack depends on BLAS/LAPACK,
  # so BLAS/LAPACK must come last.
  target_link_libraries( ${PROJECTNAME} ${EXTRA_LIB} ${BLAS_LAPACK_LIB} gomp)
elseif (BLAS_LAPACK_LIB)
  target_link_libraries( ${PROJECTNAME} ${BLAS_LAPACK_LIB} gomp)
endif()

# Install
# -------
install(TARGETS ${PROJECTNAME} RUNTIME DESTINATION bin)
