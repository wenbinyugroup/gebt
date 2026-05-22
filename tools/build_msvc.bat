@echo on
@REM Build GEBT with gfortran (cross-platform, via CMake + MinGW Makefiles).
@REM Requires gfortran on PATH. Install via:
@REM   winget install BrechtSanders.WinLibs.POSIX.UCRT
@REM For system BLAS/LAPACK/ARPACK, install MSYS2 packages:
@REM   pacman -S mingw-w64-ucrt-x86_64-openblas mingw-w64-ucrt-x86_64-arpack
@REM Then set CMAKE_PREFIX_PATH below to your MSYS2 ucrt64 lib path.

cmake -B build ^
  -DCMAKE_Fortran_COMPILER=gfortran ^
  -DCMAKE_VERBOSE_MAKEFILE=1 ^
  -DCMAKE_PREFIX_PATH=C:/msys64/ucrt64 ^
  -G "MinGW Makefiles" .

cmake --build build
