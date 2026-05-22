@echo on
@REM Build GEBT with Intel Fortran (ifx). Run from Intel oneAPI command prompt.
@REM Requires ifx on PATH.

cmake -B build_ifx ^
  -DCMAKE_Fortran_COMPILER=ifx ^
  -DCMAKE_VERBOSE_MAKEFILE=1 ^
  . -G "Visual Studio 17 2022" -T "fortran=ifx"

cmake --build build_ifx --config Release

@REM Copy Intel Fortran runtime DLLs next to gebt.exe so it can run
@REM outside the Intel oneAPI command prompt.
for /F "delims=" %%I in ('where ifx 2^>nul') do set "IFX_PATH=%%~dpI"
if defined IFX_PATH (
  for %%D in (libifcoremd.dll libifportmd.dll libmmd.dll svml_dispmd.dll) do (
    if exist "%IFX_PATH%%%D" (
      copy /Y "%IFX_PATH%%%D" "build_ifx\Release\" >nul
    ) else (
      echo WARNING: %IFX_PATH%%%D not found
    )
  )
) else (
  echo WARNING: ifx not found on PATH — skipping DLL copy
)
