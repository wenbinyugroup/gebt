# GEBT

Geometrically Exact Beam Theory solver.

Examples are provided in `examples/`.

User manuals and supporting documents are in `doc/`.

Online documentation: <https://wenbinyugroup.github.io/gebt>

## Build From Source

### Requirements

- CMake 3.10 or newer
- A Fortran compiler
  - Primary: `gfortran`
  - Fallback on Windows: Intel `ifx`
- MA28/MC19/ddep are bundled in `src/`
- BLAS/LAPACK/ARPACK
  - Preferred: system libraries
  - Fallback: bundled `src/archive/blas.f`, `src/archive/lapack.f`, `src/archive/arpack.f`

### Linux

Build with CMake:

```sh
cmake -B build .
cmake --build build
```

If system BLAS/LAPACK/ARPACK are installed in standard locations, CMake will
link them automatically. Otherwise it falls back to the bundled Fortran
sources in `src/archive/`.

### Windows with gfortran (recommended)

Install:

- `gfortran` on `PATH`
  - Example: `winget install BrechtSanders.WinLibs.POSIX.UCRT`
- MSYS2 UCRT64 OpenBLAS and ARPACK:
  - `pacman -S mingw-w64-ucrt-x86_64-openblas mingw-w64-ucrt-x86_64-arpack`

Then run:

```bat
.\tools\build_msvc.bat
```

That script configures CMake with:

- `gfortran` as the compiler
- `MinGW Makefiles` as the generator
- `CMAKE_PREFIX_PATH=C:/msys64/ucrt64` for MSYS2 system libraries

The gfortran build statically links `libgfortran` and `libgcc` on Windows to
avoid a UCRT heap conflict observed with `TRIM()` temporaries.

### Windows with Intel ifx (fallback)

Open an Intel oneAPI command prompt, go to the repository root, then run:

```bat
.\tools\build_ifx.bat
```

### Legacy Makefile

The legacy build is still available:

```sh
cd tools
make install
```
