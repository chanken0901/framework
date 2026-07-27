# ============================================================
# Makefile for current NSE/GPE common-module environment
# Windows + gfortran + Microsoft MPI
# ============================================================

FC = gfortran

# ---- Microsoft MPI ----
MPI_ROOT    = C:/Program Files (x86)/Microsoft SDKs/MPI
MPI_INC     = $(MPI_ROOT)/Include
MPI_INC_X64 = $(MPI_ROOT)/Include/x64
MPI_LIB     = $(MPI_ROOT)/Lib/x64

# ---- Directories ----
SRC_ROOT = src
OBJ_ROOT = build
MOD_DIR  = $(OBJ_ROOT)/mod
EXE      = solver.exe

# ---- Flags ----
FFLAGS_BASE = -fopenmp -fallow-invalid-boz -fallow-argument-mismatch
FFLAGS_OPT  = -O2
FFLAGS_DBG  = -O0 -g -Wall -Wextra -fcheck=all -fbacktrace
FFLAGS      = $(FFLAGS_BASE) $(FFLAGS_OPT)

INCLUDES = -I. -I$(SRC_ROOT) -I$(MOD_DIR) -I"$(MPI_INC)" -I"$(MPI_INC_X64)"
LIBS     = -L"$(MPI_LIB)" -lmsmpifec -lmsmpi -fopenmp

# ---- Source files ----
# Important: Fortran modules must be compiled in dependency order.
OBJS = \
	$(OBJ_ROOT)/common/mod_precision.o \
	$(OBJ_ROOT)/common/mod_constants.o \
	$(OBJ_ROOT)/common/mod_common_config.o \
	$(OBJ_ROOT)/common/mod_model_config.o \
	$(OBJ_ROOT)/grid/mod_grid_fvm.o \
	$(OBJ_ROOT)/field/mod_field_nse.o \
	$(OBJ_ROOT)/io/mod_input_reader.o \
	$(OBJ_ROOT)/io/mod_slf_output.o \
	$(OBJ_ROOT)/mpi/module_mpi_TW_NEC.o \
	$(OBJ_ROOT)/main/main_nse.o

.PHONY: all debug run clean dirs

all: dirs $(EXE)

# Debug build: make debug
debug:
	$(MAKE) FFLAGS="$(FFLAGS_BASE) $(FFLAGS_DBG)" all

$(EXE): $(OBJS)
	$(FC) $(OBJS) $(LIBS) -o $(EXE)

# ---- Directory creation ----
dirs:
	@if not exist "$(OBJ_ROOT)" mkdir "$(OBJ_ROOT)"
	@if not exist "$(MOD_DIR)" mkdir "$(MOD_DIR)"
	@if not exist "$(OBJ_ROOT)\common" mkdir "$(OBJ_ROOT)\common"
	@if not exist "$(OBJ_ROOT)\grid" mkdir "$(OBJ_ROOT)\grid"
	@if not exist "$(OBJ_ROOT)\field" mkdir "$(OBJ_ROOT)\field"
	@if not exist "$(OBJ_ROOT)\io" mkdir "$(OBJ_ROOT)\io"
	@if not exist "$(OBJ_ROOT)\mpi" mkdir "$(OBJ_ROOT)\mpi"
	@if not exist "$(OBJ_ROOT)\main" mkdir "$(OBJ_ROOT)\main"

# ---- Common modules ----
$(OBJ_ROOT)/common/mod_precision.o: $(SRC_ROOT)/common/mod_precision.f90 | dirs
	$(FC) -c $< $(INCLUDES) $(FFLAGS) -J$(MOD_DIR) -o $@

$(OBJ_ROOT)/common/mod_constants.o: $(SRC_ROOT)/common/mod_constants.f90 $(OBJ_ROOT)/common/mod_precision.o | dirs
	$(FC) -c $< $(INCLUDES) $(FFLAGS) -J$(MOD_DIR) -o $@

$(OBJ_ROOT)/common/mod_common_config.o: $(SRC_ROOT)/common/mod_common_config.f90 $(OBJ_ROOT)/common/mod_precision.o | dirs
	$(FC) -c $< $(INCLUDES) $(FFLAGS) -J$(MOD_DIR) -o $@

$(OBJ_ROOT)/common/mod_model_config.o: $(SRC_ROOT)/common/mod_model_config.f90 $(OBJ_ROOT)/common/mod_precision.o | dirs
	$(FC) -c $< $(INCLUDES) $(FFLAGS) -J$(MOD_DIR) -o $@

# ---- Grid modules ----
$(OBJ_ROOT)/grid/mod_grid_fvm.o: $(SRC_ROOT)/grid/mod_grid_fvm.f90 \
	$(OBJ_ROOT)/common/mod_common_config.o \
	$(OBJ_ROOT)/common/mod_precision.o | dirs
	$(FC) -c $< $(INCLUDES) $(FFLAGS) -J$(MOD_DIR) -o $@

# ---- Field modules ----
$(OBJ_ROOT)/field/mod_field_nse.o: $(SRC_ROOT)/field/mod_field_nse.f90 \
	$(OBJ_ROOT)/common/mod_precision.o \
	$(OBJ_ROOT)/common/mod_model_config.o \
	$(OBJ_ROOT)/common/mod_common_config.o | dirs
	$(FC) -c $< $(INCLUDES) $(FFLAGS) -J$(MOD_DIR) -o $@


# ---- IO module ----
$(OBJ_ROOT)/io/mod_input_reader.o: $(SRC_ROOT)/io/mod_input_reader.f90 \
	$(OBJ_ROOT)/common/mod_precision.o \
	$(OBJ_ROOT)/common/mod_common_config.o \
	$(OBJ_ROOT)/common/mod_model_config.o | dirs
	$(FC) -c $< $(INCLUDES) $(FFLAGS) -J$(MOD_DIR) -o $@

$(OBJ_ROOT)/io/mod_slf_output.o: $(SRC_ROOT)/io/mod_slf_output.f90 \
	$(OBJ_ROOT)/common/mod_precision.o \
	$(OBJ_ROOT)/common/mod_common_config.o \
	$(OBJ_ROOT)/common/mod_model_config.o \
	$(OBJ_ROOT)/mpi/module_mpi_TW_NEC.o | dirs
	$(FC) -c $< $(INCLUDES) $(FFLAGS) -J$(MOD_DIR) -o $@

# ---- MPI module ----
$(OBJ_ROOT)/mpi/module_mpi_TW_NEC.o: $(SRC_ROOT)/mpi/module_mpi_TW_NEC.f90 | dirs
	$(FC) -c $< $(INCLUDES) $(FFLAGS) -J$(MOD_DIR) -o $@


# ---- Main program ----
$(OBJ_ROOT)/main/main_nse.o: $(SRC_ROOT)/main/main_nse.f90 \
	$(OBJ_ROOT)/common/mod_precision.o \
	$(OBJ_ROOT)/common/mod_constants.o \
	$(OBJ_ROOT)/common/mod_common_config.o \
	$(OBJ_ROOT)/common/mod_model_config.o \
	$(OBJ_ROOT)/io/mod_input_reader.o \
	$(OBJ_ROOT)/io/mod_slf_output.o \
	$(OBJ_ROOT)/mpi/module_mpi_TW_NEC.o | dirs
	$(FC) -c $< $(INCLUDES) $(FFLAGS) -J$(MOD_DIR) -o $@

# ---- Execution ----
# Usage: make run NP=4
NP ?= 4
OMP ?= 8
run: all
	mpiexec -n $(NP) $(EXE)

# ---- Clean ----
clean:
	-del /Q $(EXE) 2>nul
	-del /Q *.mod 2>nul
	-if exist $(OBJ_ROOT) rmdir /S /Q $(OBJ_ROOT)
