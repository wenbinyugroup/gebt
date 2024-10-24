@echo on

mkdir build_mingw_64
cd build_mingw_64
cmake ^
  -DCMAKE_C_COMPILER=x86_64-w64-mingw32-gcc ^
  -DCMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++ ^
  -DCMAKE_Fortran_COMPILER=x86_64-w64-mingw32-gfortran ^
  -DCMAKE_PREFIX_PATH="C:/cygwin64/usr/x86_64-w64-mingw32/sys-root/mingw" ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_INSTALL_PREFIX=.. ^
  -DCMAKE_VERBOSE_MAKEFILE=1 ^
  .. -G"Unix Makefiles"

make install

cd ..

