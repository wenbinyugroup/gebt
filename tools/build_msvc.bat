@echo on

@REM mkdir build_msvc
@REM cd build_msvc

cmake -B build_msvc ^
  -DCMAKE_Fortran_COMPILER=ifx ^
  -DCMAKE_INSTALL_PREFIX=.. ^
  -DCMAKE_VERBOSE_MAKEFILE=1 ^
  . -G "Visual Studio 17 2022" -T "fortran=ifx"

cmake --build build_msvc --config Release

@REM cd ..

