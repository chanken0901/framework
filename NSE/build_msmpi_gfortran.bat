@echo off
setlocal enabledelayedexpansion

rem ============================================================
rem build_msmpi_gfortran.bat
rem
rem MS-MPI + gfortran 用ビルドスクリプト
rem
rem 対象:
rem   - module_mpi_TW_NEC.f90
rem   - 3d_mpi_openmp_ver2_TGV_KEEP.f90
rem
rem 使い方:
rem   1. このbatを Fortran ファイルと同じフォルダに置く
rem   2. cmd でそのフォルダへ移動
rem   3. build_msmpi_gfortran.bat
rem
rem 実行:
rem   mpiexec -n 4 solver.exe
rem ============================================================

rem ---------- User settings ----------
set FC=gfortran

set MPI_ROOT=C:\Program Files (x86)\Microsoft SDKs\MPI
set MPI_INC=%MPI_ROOT%\Include
set MPI_INC_X64=%MPI_ROOT%\Include\x64
set MPI_LIB=%MPI_ROOT%\Lib\x64

set MPI_MODULE_SRC=module_mpi_TW_NEC.f90
set MAIN_SRC=3d_mpi_openmp_ver2_TGV_KEEP.f90
set EXE=solver.exe

rem ---------- Compiler / linker flags ----------
rem -fopenmp:
rem   OpenMP を有効化
rem
rem -fallow-invalid-boz:
rem   MS-MPI の mpif.h 内の BOZ literal 警告/エラー対策
rem
rem -fallow-argument-mismatch:
rem   mpif.h を使う古いMPIコードに対する引数型不一致対策
rem
set FFLAGS=-O0 -fopenmp -fallow-invalid-boz -fallow-argument-mismatch
set INCLUDES=-I. -I"%MPI_INC%" -I"%MPI_INC_X64%"
set LIBS=-L"%MPI_LIB%" -lmsmpifec -lmsmpi -fopenmp

echo.
echo ============================================================
echo  MS-MPI + gfortran build
echo ============================================================
echo FC      = %FC%
echo MPI_ROOT= %MPI_ROOT%
echo MODULE  = %MPI_MODULE_SRC%
echo MAIN    = %MAIN_SRC%
echo EXE     = %EXE%
echo.

rem ---------- Check files ----------
if not exist "%MPI_MODULE_SRC%" (
    echo [ERROR] %MPI_MODULE_SRC% not found.
    goto :error
)

if not exist "%MAIN_SRC%" (
    echo [ERROR] %MAIN_SRC% not found.
    goto :error
)

if not exist "%MPI_INC%\mpif.h" (
    echo [ERROR] mpif.h not found:
    echo         %MPI_INC%\mpif.h
    goto :error
)

if not exist "%MPI_LIB%\msmpi.lib" (
    echo [ERROR] msmpi.lib not found:
    echo         %MPI_LIB%\msmpi.lib
    goto :error
)

if not exist "%MPI_LIB%\msmpifec.lib" (
    echo [ERROR] msmpifec.lib not found:
    echo         %MPI_LIB%\msmpifec.lib
    echo.
    echo Check:
    echo   dir "%MPI_LIB%"
    goto :error
)

rem ---------- Clean old objects ----------
echo [1/4] Cleaning old build files...
del /Q *.o *.mod "%EXE%" 2>nul

rem ---------- Compile MPI module ----------
echo.
echo [2/4] Compiling MPI module...
%FC% -c "%MPI_MODULE_SRC%" %INCLUDES% %FFLAGS%
if errorlevel 1 (
    echo.
    echo [ERROR] Failed to compile %MPI_MODULE_SRC%.
    goto :error
)

if not exist "module_mpi.mod" (
    echo.
    echo [ERROR] module_mpi.mod was not generated.
    echo         Check module name in %MPI_MODULE_SRC%.
    goto :error
)

rem ---------- Compile main ----------
echo.
echo [3/4] Compiling main program...
%FC% -c "%MAIN_SRC%" %INCLUDES% %FFLAGS%
if errorlevel 1 (
    echo.
    echo [ERROR] Failed to compile %MAIN_SRC%.
    goto :error
)

rem ---------- Link ----------
echo.
echo [4/4] Linking...
%FC% module_mpi_TW_NEC.o 3d_mpi_openmp_ver2_TGV_KEEP.o %LIBS% -o "%EXE%"
if errorlevel 1 (
    echo.
    echo [ERROR] Failed to link %EXE%.
    goto :error
)

echo.
echo ============================================================
echo [OK] Build completed.
echo ============================================================
echo Output:
echo   %EXE%
echo.
echo Example run:
echo   mpiexec -n 4 %EXE%
echo.
goto :end

:error
echo.
echo ============================================================
echo [FAILED] Build failed.
echo ============================================================
echo.
echo Useful checks:
echo   gfortran --version
echo   mpiexec
echo   dir "%MPI_INC%"
echo   dir "%MPI_LIB%"
echo.

:end
endlocal
