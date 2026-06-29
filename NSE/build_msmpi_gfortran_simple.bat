@echo off
echo ============================================
echo Compile module_mpi_TW_NEC.f90
echo ============================================

gfortran -c module_mpi_TW_NEC.f90 -I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include" -I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include\x64" -fopenmp -fallow-invalid-boz -fallow-argument-mismatch
if errorlevel 1 goto :error

echo.
echo ============================================
echo Compile main
echo ============================================

gfortran -c 3d_mpi_openmp_ver2_TGV_KEEP.f90 -I. -I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include" -I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include\x64" -fopenmp -fallow-invalid-boz -fallow-argument-mismatch
if errorlevel 1 goto :error

echo.
echo ============================================
echo Link
echo ============================================

gfortran module_mpi_TW_NEC.o 3d_mpi_openmp_ver2_TGV_KEEP.o -L"C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" -lmsmpifec -lmsmpi -fopenmp -o solver.exe
if errorlevel 1 goto :error

echo.
echo Build completed successfully.
echo.
echo Run:
echo mpiexec -n 4 solver.exe
pause
exit /b

:error
echo.
echo Build failed.
pause
