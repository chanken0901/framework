#include <cuda_runtime.h>
#include <cufftMp.h>
#include <cuComplex.h>
#include <mpi.h>

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <utility>

#define GP3D_CUFFTMP_API extern "C"

namespace {

constexpr int kThreads = 256;
constexpr double kPi = 3.141592653589793238462643383279502884;
char g_last_error[1024] = "";

struct GpuContext {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  int rank = 0;
  int nprocs = 1;
  int local_rank = 0;
  int local_size = 1;
  int device = 0;
  int z_start = 0;
  int local_nz = 0;
  int y_start = 0;
  int local_ny = 0;
  size_t n_global = 0;
  size_t n_real = 0;
  size_t n_spectral = 0;
  size_t workspace_bytes = 0;
  bool has_argle = false;
  bool plan_created = false;
  MPI_Comm comm = MPI_COMM_NULL;
  MPI_Comm local_comm = MPI_COMM_NULL;
  cufftHandle plan{};
  cudaLibXtDesc* psi_desc = nullptr;
  cudaLibXtDesc* work_desc = nullptr;
  cudaLibXtDesc* old_psi_desc = nullptr;
  cudaLibXtDesc* grad_x_desc = nullptr;
  cudaLibXtDesc* grad_y_desc = nullptr;
  cudaLibXtDesc* rhs_desc = nullptr;
  double* potential = nullptr;
  double* kx = nullptr;
  double* ky = nullptr;
  double* kz = nullptr;
  double* reductions = nullptr;
};

void set_error(const char* message) {
  std::snprintf(g_last_error, sizeof(g_last_error), "%s", message);
}

void set_error(const char* operation, cudaError_t error) {
  std::snprintf(g_last_error, sizeof(g_last_error), "%s: %s", operation,
                cudaGetErrorString(error));
}

void set_error(const char* operation, cufftResult error) {
  std::snprintf(g_last_error, sizeof(g_last_error), "%s: cuFFTMp error %d",
                operation, static_cast<int>(error));
}

void set_error(const char* operation, int mpi_error) {
  char mpi_message[MPI_MAX_ERROR_STRING] = "";
  int length = 0;
  MPI_Error_string(mpi_error, mpi_message, &length);
  std::snprintf(g_last_error, sizeof(g_last_error), "%s: MPI error %.*s",
                operation, length, mpi_message);
}

bool cuda_ok(cudaError_t error, const char* operation) {
  if (error == cudaSuccess) return true;
  set_error(operation, error);
  return false;
}

bool cufft_ok(cufftResult error, const char* operation) {
  if (error == CUFFT_SUCCESS) return true;
  set_error(operation, error);
  return false;
}

bool mpi_ok(int error, const char* operation) {
  if (error == MPI_SUCCESS) return true;
  set_error(operation, error);
  return false;
}

bool kernel_ok(const char* operation) {
  return cuda_ok(cudaGetLastError(), operation);
}

int block_count(size_t n) {
  return static_cast<int>((n + static_cast<size_t>(kThreads) - 1) /
                          static_cast<size_t>(kThreads));
}

void block_range(int n, int rank, int nprocs, int* start, int* count) {
  const int base = n / nprocs;
  const int rest = n % nprocs;
  if (rank < rest) {
    *count = base + 1;
    *start = rank * (*count);
  } else {
    *count = base;
    *start = rest * (base + 1) + (rank - rest) * base;
  }
}

cufftDoubleComplex* descriptor_data(cudaLibXtDesc* desc) {
  if (desc == nullptr || desc->descriptor == nullptr) return nullptr;
  return static_cast<cufftDoubleComplex*>(desc->descriptor->data[0]);
}

bool allocate_fft_buffer(GpuContext* ctx, cudaLibXtDesc** desc,
                         const char* operation) {
  if (!cufft_ok(cufftXtMalloc(ctx->plan, desc, CUFFT_XT_FORMAT_INPLACE),
                operation)) {
    return false;
  }
  if (descriptor_data(*desc) == nullptr) {
    set_error("cuFFTMp returned a descriptor without local GPU memory");
    return false;
  }
  return true;
}

void free_fft_buffer(cudaLibXtDesc*& desc) {
  if (desc != nullptr) cufftXtFree(desc);
  desc = nullptr;
}

void destroy_context(GpuContext* ctx) {
  if (ctx == nullptr) return;
  cudaSetDevice(ctx->device);
  cudaDeviceSynchronize();

  // Descriptor memory must be released while the cuFFTMp plan is alive.
  free_fft_buffer(ctx->rhs_desc);
  free_fft_buffer(ctx->grad_y_desc);
  free_fft_buffer(ctx->grad_x_desc);
  free_fft_buffer(ctx->old_psi_desc);
  free_fft_buffer(ctx->work_desc);
  free_fft_buffer(ctx->psi_desc);
  if (ctx->plan_created) cufftDestroy(ctx->plan);

  cudaFree(ctx->potential);
  cudaFree(ctx->kx);
  cudaFree(ctx->ky);
  cudaFree(ctx->kz);
  cudaFree(ctx->reductions);

  int finalized = 0;
  MPI_Finalized(&finalized);
  if (!finalized) {
    if (ctx->local_comm != MPI_COMM_NULL) MPI_Comm_free(&ctx->local_comm);
    if (ctx->comm != MPI_COMM_NULL) MPI_Comm_free(&ctx->comm);
  }
  ctx->~GpuContext();
  std::free(ctx);
}

bool allreduce_values(GpuContext* ctx, const double* local, double* global,
                      int count, MPI_Op operation, const char* label) {
  return mpi_ok(MPI_Allreduce(local, global, count, MPI_DOUBLE, operation,
                              ctx->comm),
                label);
}

bool synchronize_for_timing(GpuContext* ctx, double* timestamp) {
  if (!cuda_ok(cudaDeviceSynchronize(), "synchronize timed CUDA work")) {
    return false;
  }
  *timestamp = MPI_Wtime();
  return true;
}

__device__ cufftDoubleComplex complex_scale(cufftDoubleComplex value,
                                            double scale) {
  return make_cuDoubleComplex(value.x * scale, value.y * scale);
}

__device__ cufftDoubleComplex complex_multiply(cufftDoubleComplex a,
                                               cufftDoubleComplex b) {
  return make_cuDoubleComplex(a.x * b.x - a.y * b.y,
                             a.x * b.y + a.y * b.x);
}

__device__ double complex_abs2(cufftDoubleComplex value) {
  return value.x * value.x + value.y * value.y;
}

__device__ double atomic_add_double(double* address, double value) {
#if __CUDA_ARCH__ >= 600
  return atomicAdd(address, value);
#else
  auto* bits = reinterpret_cast<unsigned long long*>(address);
  unsigned long long old = *bits;
  unsigned long long assumed;
  do {
    assumed = old;
    old = atomicCAS(bits, assumed,
                    __double_as_longlong(value + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
#endif
}

__device__ void atomic_min_positive(double* address, double value) {
  atomicMin(reinterpret_cast<unsigned long long*>(address),
            __double_as_longlong(value));
}

__device__ void atomic_max_positive(double* address, double value) {
  atomicMax(reinterpret_cast<unsigned long long*>(address),
            __double_as_longlong(value));
}

__global__ void local_half_step_kernel(cufftDoubleComplex* psi,
                                       const double* potential, size_t n,
                                       double tau, double g,
                                       int imaginary_time) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  if (index >= n) return;
  const cufftDoubleComplex value = psi[index];
  const double energy = potential[index] + g * complex_abs2(value);
  cufftDoubleComplex factor;
  if (imaginary_time != 0) {
    factor = make_cuDoubleComplex(exp(-tau * energy), 0.0);
  } else {
    const double phase = -tau * energy;
    factor = make_cuDoubleComplex(cos(phase), sin(phase));
  }
  psi[index] = complex_multiply(value, factor);
}

// cuFFTMp's forward transform changes [local_z][ny][nx] into
// [nz][local_y][nx].  This kernel works directly on that shuffled layout.
__global__ void kinetic_shuffled_kernel(
    cufftDoubleComplex* spectral, const double* kx, const double* ky,
    const double* kz, int nx, int local_ny, int y_start, int nz, size_t n,
    double tau, double kinetic_coefficient, int imaginary_time) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  if (index >= n) return;
  const int i = static_cast<int>(index % nx);
  const size_t row = index / nx;
  const int j = static_cast<int>(row % local_ny) + y_start;
  const int k = static_cast<int>(row / local_ny);
  const double k2 = kx[i] * kx[i] + ky[j] * ky[j] + kz[k] * kz[k];
  const double energy = kinetic_coefficient * k2;
  cufftDoubleComplex factor;
  if (imaginary_time != 0) {
    factor = make_cuDoubleComplex(exp(-tau * energy), 0.0);
  } else {
    const double phase = -tau * energy;
    factor = make_cuDoubleComplex(cos(phase), sin(phase));
  }
  spectral[index] = complex_multiply(spectral[index], factor);
}

__global__ void scale_kernel(cufftDoubleComplex* values, size_t n,
                             double scale) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  if (index < n) values[index] = complex_scale(values[index], scale);
}

__global__ void norm_kernel(const cufftDoubleComplex* psi, size_t n,
                            double* sum) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  if (index < n) atomic_add_double(sum, complex_abs2(psi[index]));
}

__global__ void derivative_shuffled_kernel(
    const cufftDoubleComplex* spectral, cufftDoubleComplex* derivative,
    const double* kx, const double* ky, const double* kz, int nx,
    int local_ny, int y_start, size_t n, int axis) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  if (index >= n) return;
  const int i = static_cast<int>(index % nx);
  const size_t row = index / nx;
  const int j = static_cast<int>(row % local_ny) + y_start;
  const int k = static_cast<int>(row / local_ny);
  const double wave_number = axis == 0 ? kx[i] : (axis == 1 ? ky[j] : kz[k]);
  const cufftDoubleComplex value = spectral[index];
  derivative[index] = make_cuDoubleComplex(-wave_number * value.y,
                                           wave_number * value.x);
}

__global__ void argle_rhs_natural_kernel(
    const cufftDoubleComplex* old_psi, const cufftDoubleComplex* grad_x,
    const cufftDoubleComplex* grad_y, cufftDoubleComplex* rhs,
    const double* potential, int nx, int ny, int local_nz, int z_start,
    int nz, size_t n, double reaction_constant,
    double nonlinear_coefficient, double potential_scale, double alpha,
    double velocity_amplitude) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  if (index >= n) return;
  const int i = static_cast<int>(index % nx);
  const size_t plane = index / nx;
  const int j = static_cast<int>(plane % ny);
  const int k = static_cast<int>(plane / ny) + z_start;

  const double x = 2.0 * kPi * (static_cast<double>(i) + 0.5) / nx - kPi;
  const double y = 2.0 * kPi * (static_cast<double>(j) + 0.5) / ny - kPi;
  const double z = 2.0 * kPi * (static_cast<double>(k) + 0.5) / nz - kPi;
  const double vx = velocity_amplitude * sin(x) * cos(y) * cos(z);
  const double vy = -velocity_amplitude * cos(x) * sin(y) * cos(z);
  const double velocity2 = vx * vx + vy * vy;

  const cufftDoubleComplex value = old_psi[index];
  const double reaction = reaction_constant - potential_scale * potential[index] -
                          nonlinear_coefficient * complex_abs2(value) -
                          velocity2 / (4.0 * alpha);
  const cufftDoubleComplex advection = make_cuDoubleComplex(
      vx * grad_x[index].x + vy * grad_y[index].x,
      vx * grad_x[index].y + vy * grad_y[index].y);
  rhs[index] = make_cuDoubleComplex(reaction * value.x + advection.y,
                                    reaction * value.y - advection.x);
}

__global__ void argle_update_shuffled_kernel(
    cufftDoubleComplex* spectral, const cufftDoubleComplex* rhs_spectral,
    const double* kx, const double* ky, const double* kz, int nx,
    int local_ny, int y_start, size_t n, double alpha, double dtau) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  if (index >= n) return;
  const int i = static_cast<int>(index % nx);
  const size_t row = index / nx;
  const int j = static_cast<int>(row % local_ny) + y_start;
  const int k = static_cast<int>(row / local_ny);
  const double k2 = kx[i] * kx[i] + ky[j] * ky[j] + kz[k] * kz[k];
  const double implicit_term = 0.5 * dtau * alpha * k2;
  const double denominator = 1.0 + implicit_term;
  const double explicit_factor = 1.0 - implicit_term;
  spectral[index] = make_cuDoubleComplex(
      (explicit_factor * spectral[index].x + dtau * rhs_spectral[index].x) /
          denominator,
      (explicit_factor * spectral[index].y + dtau * rhs_spectral[index].y) /
          denominator);
}

__global__ void argle_metrics_kernel(const cufftDoubleComplex* psi,
                                     const cufftDoubleComplex* old_psi,
                                     size_t n, double* reductions) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  if (index >= n) return;
  const double density = complex_abs2(psi[index]);
  const double dr = psi[index].x - old_psi[index].x;
  const double di = psi[index].y - old_psi[index].y;
  atomic_add_double(&reductions[0], density);
  atomic_min_positive(&reductions[1], density);
  atomic_max_positive(&reductions[2], sqrt(dr * dr + di * di));
}

__global__ void local_diagnostics_kernel(const cufftDoubleComplex* psi,
                                         const double* potential, size_t n,
                                         double g, double* reductions) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  if (index >= n) return;
  const double density = complex_abs2(psi[index]);
  atomic_add_double(&reductions[0], density);
  atomic_add_double(&reductions[1],
                    potential[index] * density + 0.5 * g * density * density);
}

__global__ void spectral_diagnostics_shuffled_kernel(
    const cufftDoubleComplex* spectral, const double* kx, const double* ky,
    const double* kz, int nx, int local_ny, int y_start, size_t n,
    double* reductions) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  if (index >= n) return;
  const int i = static_cast<int>(index % nx);
  const size_t row = index / nx;
  const int j = static_cast<int>(row % local_ny) + y_start;
  const int k = static_cast<int>(row / local_ny);
  const double k2 = kx[i] * kx[i] + ky[j] * ky[j] + kz[k] * kz[k];
  atomic_add_double(&reductions[2], k2 * complex_abs2(spectral[index]));
}

bool forward_fft(GpuContext* ctx, cufftDoubleComplex* data,
                 const char* operation) {
  return cufft_ok(cufftExecZ2Z(ctx->plan, data, data, CUFFT_FORWARD), operation);
}

bool inverse_fft(GpuContext* ctx, cufftDoubleComplex* data,
                 const char* operation) {
  return cufft_ok(cufftExecZ2Z(ctx->plan, data, data, CUFFT_INVERSE), operation);
}

bool normalize_wavefunction(GpuContext* ctx, double target_norm,
                            double volume_element) {
  if (!cuda_ok(cudaMemset(ctx->reductions, 0, sizeof(double)),
               "clear local norm reduction")) {
    return false;
  }
  auto* psi = descriptor_data(ctx->psi_desc);
  norm_kernel<<<block_count(ctx->n_real), kThreads>>>(psi, ctx->n_real,
                                                      ctx->reductions);
  if (!kernel_ok("local norm kernel")) return false;
  double local_sum = 0.0;
  if (!cuda_ok(cudaMemcpy(&local_sum, ctx->reductions, sizeof(double),
                          cudaMemcpyDeviceToHost),
               "copy local norm reduction")) {
    return false;
  }
  double global_sum = 0.0;
  if (!allreduce_values(ctx, &local_sum, &global_sum, 1, MPI_SUM,
                        "allreduce wave-function norm")) {
    return false;
  }
  const double current_norm = global_sum * volume_element;
  if (!(current_norm > 0.0)) {
    set_error("cannot normalize a zero distributed wave function");
    return false;
  }
  const double scale = sqrt(target_norm / current_norm);
  scale_kernel<<<block_count(ctx->n_real), kThreads>>>(psi, ctx->n_real, scale);
  return kernel_ok("distributed normalization scale kernel");
}

void print_topology(GpuContext* ctx) {
  cudaDeviceProp properties{};
  int cufft_version = 0;
  int cuda_runtime_version = 0;
  int cuda_driver_version = 0;
  cudaGetDeviceProperties(&properties, ctx->device);
  cufftGetVersion(&cufft_version);
  cudaRuntimeGetVersion(&cuda_runtime_version);
  cudaDriverGetVersion(&cuda_driver_version);
  if (ctx->rank == 0) {
#if defined(GP3D_CUFFTMP_LEGACY_API)
    const char* api_name = "legacy attach-communicator API";
#else
    const char* api_name = "direct communicator plan API";
#endif
    std::printf("# GPU backend: MPI + cuFFTMp (%s)\n", api_name);
    std::printf("# cuFFTMp ranks=%d global_grid=%dx%dx%d workspace_bytes=%zu\n",
                ctx->nprocs, ctx->nx, ctx->ny, ctx->nz,
                ctx->workspace_bytes);
    std::printf("# library_versions cufftMp=%d cuda_runtime=%d cuda_driver=%d\n",
                cufft_version, cuda_runtime_version, cuda_driver_version);
  }
  for (int owner = 0; owner < ctx->nprocs; ++owner) {
    MPI_Barrier(ctx->comm);
    if (owner == ctx->rank) {
      std::printf(
          "# cufftmp rank=%d local_rank=%d device=%d name=\"%s\" "
          "real_z=[%d,%d] spectral_y=[%d,%d]\n",
          ctx->rank, ctx->local_rank, ctx->device, properties.name,
          ctx->z_start + 1, ctx->z_start + ctx->local_nz,
          ctx->y_start + 1, ctx->y_start + ctx->local_ny);
      std::fflush(stdout);
    }
  }
  MPI_Barrier(ctx->comm);
}

}  // namespace

GP3D_CUFFTMP_API int gp3d_cufftmp_create(
    void** handle, int nx, int ny, int nz, int expected_local_nz,
    int expected_k_start, int need_argle, int fortran_comm) {
  g_last_error[0] = '\0';
  if (handle == nullptr || nx <= 0 || ny <= 0 || nz <= 0) {
    set_error("invalid cuFFTMp context dimensions");
    return 1;
  }
  *handle = nullptr;

  int initialized = 0;
  MPI_Initialized(&initialized);
  if (!initialized) {
    set_error("MPI must be initialized before creating the cuFFTMp context");
    return 1;
  }

  void* storage = std::malloc(sizeof(GpuContext));
  if (storage == nullptr) {
    set_error("failed to allocate cuFFTMp context");
    return 1;
  }
  auto* ctx = ::new (storage) GpuContext();
  ctx->nx = nx;
  ctx->ny = ny;
  ctx->nz = nz;
  ctx->has_argle = need_argle != 0;

  const MPI_Comm source_comm =
      MPI_Comm_f2c(static_cast<MPI_Fint>(fortran_comm));
  if (!mpi_ok(MPI_Comm_dup(source_comm, &ctx->comm),
              "duplicate cuFFTMp communicator")) {
    destroy_context(ctx);
    return 1;
  }
  if (!mpi_ok(MPI_Comm_rank(ctx->comm, &ctx->rank), "query MPI rank") ||
      !mpi_ok(MPI_Comm_size(ctx->comm, &ctx->nprocs), "query MPI size")) {
    destroy_context(ctx);
    return 1;
  }
  if (ctx->nprocs > nz || ctx->nprocs > ny) {
    set_error("cuFFTMp built-in slabs require both nz and ny to be at least the MPI process count");
    destroy_context(ctx);
    return 1;
  }

  if (!mpi_ok(MPI_Comm_split_type(ctx->comm, MPI_COMM_TYPE_SHARED, ctx->rank,
                                  MPI_INFO_NULL, &ctx->local_comm),
              "create node-local MPI communicator") ||
      !mpi_ok(MPI_Comm_rank(ctx->local_comm, &ctx->local_rank),
              "query node-local MPI rank") ||
      !mpi_ok(MPI_Comm_size(ctx->local_comm, &ctx->local_size),
              "query node-local MPI size")) {
    destroy_context(ctx);
    return 1;
  }

  int device_count = 0;
  if (!cuda_ok(cudaGetDeviceCount(&device_count), "query visible CUDA devices")) {
    destroy_context(ctx);
    return 1;
  }
  if (device_count <= 0 || ctx->local_size > device_count) {
    std::snprintf(g_last_error, sizeof(g_last_error),
                  "node has %d MPI ranks but only %d visible CUDA devices; "
                  "use one rank per GPU and identical CUDA_VISIBLE_DEVICES on the node",
                  ctx->local_size, device_count);
    destroy_context(ctx);
    return 1;
  }
  ctx->device = ctx->local_rank;
  if (!cuda_ok(cudaSetDevice(ctx->device), "select rank-local CUDA device")) {
    destroy_context(ctx);
    return 1;
  }

  block_range(nz, ctx->rank, ctx->nprocs, &ctx->z_start, &ctx->local_nz);
  block_range(ny, ctx->rank, ctx->nprocs, &ctx->y_start, &ctx->local_ny);
  if (ctx->local_nz != expected_local_nz ||
      ctx->z_start + 1 != expected_k_start) {
    set_error("Fortran grid decomposition does not match cuFFTMp's natural z-slab layout");
    destroy_context(ctx);
    return 1;
  }
  ctx->n_global = static_cast<size_t>(nx) * ny * nz;
  ctx->n_real = static_cast<size_t>(nx) * ny * ctx->local_nz;
  ctx->n_spectral = static_cast<size_t>(nx) * ctx->local_ny * nz;

  if (!cufft_ok(cufftCreate(&ctx->plan), "create cuFFTMp plan handle")) {
    destroy_context(ctx);
    return 1;
  }
  ctx->plan_created = true;
#if defined(GP3D_CUFFTMP_LEGACY_API)
  if (!cufft_ok(cufftMpAttachComm(ctx->plan, CUFFT_COMM_MPI, &ctx->comm),
                "attach MPI communicator to cuFFTMp plan") ||
      !cufft_ok(cufftMakePlan3d(ctx->plan, nz, ny, nx, CUFFT_Z2Z,
                               &ctx->workspace_bytes),
                "create legacy distributed Z2Z plan")) {
    destroy_context(ctx);
    return 1;
  }
#else
  if (!cufft_ok(cufftMpMakePlan3d(ctx->plan, nz, ny, nx, CUFFT_Z2Z,
                                 &ctx->comm, CUFFT_COMM_MPI,
                                 &ctx->workspace_bytes),
                "create distributed Z2Z plan")) {
    destroy_context(ctx);
    return 1;
  }
#endif
  if (!cufft_ok(cufftXtSetSubformatDefault(
                    ctx->plan, CUFFT_XT_FORMAT_INPLACE,
                    CUFFT_XT_FORMAT_INPLACE_SHUFFLED),
                "set natural and shuffled cuFFTMp layouts")) {
    destroy_context(ctx);
    return 1;
  }

  if (!allocate_fft_buffer(ctx, &ctx->psi_desc,
                           "allocate distributed psi buffer") ||
      !allocate_fft_buffer(ctx, &ctx->work_desc,
                           "allocate distributed FFT work buffer")) {
    destroy_context(ctx);
    return 1;
  }
  if (ctx->has_argle &&
      (!allocate_fft_buffer(ctx, &ctx->old_psi_desc,
                            "allocate distributed ARGLE old-psi buffer") ||
       !allocate_fft_buffer(ctx, &ctx->grad_x_desc,
                            "allocate distributed ARGLE x-gradient buffer") ||
       !allocate_fft_buffer(ctx, &ctx->grad_y_desc,
                            "allocate distributed ARGLE y-gradient buffer") ||
       !allocate_fft_buffer(ctx, &ctx->rhs_desc,
                            "allocate distributed ARGLE rhs buffer"))) {
    destroy_context(ctx);
    return 1;
  }

#define GP3D_CUDA_ALLOC(member, bytes, label)                                  \
  do {                                                                         \
    if (!cuda_ok(cudaMalloc(reinterpret_cast<void**>(&(member)), (bytes)),      \
                 (label))) {                                                   \
      destroy_context(ctx);                                                    \
      return 1;                                                                \
    }                                                                          \
  } while (false)

  GP3D_CUDA_ALLOC(ctx->potential, ctx->n_real * sizeof(double),
                  "allocate local GPU potential");
  GP3D_CUDA_ALLOC(ctx->kx, static_cast<size_t>(nx) * sizeof(double),
                  "allocate GPU kx");
  GP3D_CUDA_ALLOC(ctx->ky, static_cast<size_t>(ny) * sizeof(double),
                  "allocate GPU ky");
  GP3D_CUDA_ALLOC(ctx->kz, static_cast<size_t>(nz) * sizeof(double),
                  "allocate GPU kz");
  GP3D_CUDA_ALLOC(ctx->reductions, 3 * sizeof(double),
                  "allocate GPU reduction buffer");
#undef GP3D_CUDA_ALLOC

  if (!cuda_ok(cudaDeviceSynchronize(), "finish cuFFTMp initialization")) {
    destroy_context(ctx);
    return 1;
  }
  print_topology(ctx);
  *handle = ctx;
  return 0;
}

GP3D_CUFFTMP_API int gp3d_cufftmp_upload(
    void* handle, const cufftDoubleComplex* psi, const double* potential,
    const double* kx, const double* ky, const double* kz) {
  auto* ctx = static_cast<GpuContext*>(handle);
  if (ctx == nullptr) {
    set_error("cuFFTMp context is null during upload");
    return 1;
  }
  if (!cuda_ok(cudaSetDevice(ctx->device), "select CUDA device for upload")) {
    return 1;
  }
  if (!cuda_ok(cudaMemcpy(descriptor_data(ctx->psi_desc), psi,
                          ctx->n_real * sizeof(cufftDoubleComplex),
                          cudaMemcpyHostToDevice),
               "upload local psi slab") ||
      !cuda_ok(cudaMemcpy(ctx->potential, potential,
                          ctx->n_real * sizeof(double),
                          cudaMemcpyHostToDevice),
               "upload local potential slab") ||
      !cuda_ok(cudaMemcpy(ctx->kx, kx,
                          static_cast<size_t>(ctx->nx) * sizeof(double),
                          cudaMemcpyHostToDevice),
               "upload kx") ||
      !cuda_ok(cudaMemcpy(ctx->ky, ky,
                          static_cast<size_t>(ctx->ny) * sizeof(double),
                          cudaMemcpyHostToDevice),
               "upload ky") ||
      !cuda_ok(cudaMemcpy(ctx->kz, kz,
                          static_cast<size_t>(ctx->nz) * sizeof(double),
                          cudaMemcpyHostToDevice),
               "upload kz")) {
    return 1;
  }
  return 0;
}

GP3D_CUFFTMP_API int gp3d_cufftmp_download(void* handle,
                                           cufftDoubleComplex* psi) {
  auto* ctx = static_cast<GpuContext*>(handle);
  if (ctx == nullptr) {
    set_error("cuFFTMp context is null during download");
    return 1;
  }
  if (!cuda_ok(cudaSetDevice(ctx->device), "select CUDA device for download")) {
    return 1;
  }
  return cuda_ok(cudaMemcpy(psi, descriptor_data(ctx->psi_desc),
                            ctx->n_real * sizeof(cufftDoubleComplex),
                            cudaMemcpyDeviceToHost),
                 "download local psi slab")
             ? 0
             : 1;
}

GP3D_CUFFTMP_API int gp3d_cufftmp_step(
    void* handle, double dt, double hbar, double mass, double g,
    double target_norm, double volume_element, int imaginary_time,
    int measure_timing, double* nonlinear_seconds, double* fft_seconds,
    double* kinetic_seconds, double* other_seconds, double* total_seconds) {
  auto* ctx = static_cast<GpuContext*>(handle);
  if (ctx == nullptr || hbar <= 0.0 || mass <= 0.0) {
    set_error("invalid distributed CUDA split-step context or parameters");
    return 1;
  }
  if (!cuda_ok(cudaSetDevice(ctx->device), "select CUDA device for split step")) {
    return 1;
  }
  *nonlinear_seconds = 0.0;
  *fft_seconds = 0.0;
  *kinetic_seconds = 0.0;
  *other_seconds = 0.0;
  *total_seconds = 0.0;

  double stamps[8]{};
  if (measure_timing != 0) {
    if (!mpi_ok(MPI_Barrier(ctx->comm), "align timed split step") ||
        !synchronize_for_timing(ctx, &stamps[0])) {
      return 1;
    }
  }

  auto* psi = descriptor_data(ctx->psi_desc);
  const double local_tau = 0.5 * dt / hbar;
  local_half_step_kernel<<<block_count(ctx->n_real), kThreads>>>(
      psi, ctx->potential, ctx->n_real, local_tau, g, imaginary_time);
  if (!kernel_ok("first local half-step kernel")) return 1;
  if (measure_timing != 0 && !synchronize_for_timing(ctx, &stamps[1])) return 1;

  if (!forward_fft(ctx, psi, "forward distributed Z2Z transform")) return 1;
  if (measure_timing != 0 && !synchronize_for_timing(ctx, &stamps[2])) return 1;

  kinetic_shuffled_kernel<<<block_count(ctx->n_spectral), kThreads>>>(
      psi, ctx->kx, ctx->ky, ctx->kz, ctx->nx, ctx->local_ny,
      ctx->y_start, ctx->nz, ctx->n_spectral, dt / hbar,
      0.5 * hbar * hbar / mass, imaginary_time);
  if (!kernel_ok("shuffled spectral kinetic kernel")) return 1;
  if (measure_timing != 0 && !synchronize_for_timing(ctx, &stamps[3])) return 1;

  if (!inverse_fft(ctx, psi, "inverse distributed Z2Z transform")) return 1;
  if (measure_timing != 0 && !synchronize_for_timing(ctx, &stamps[4])) return 1;

  scale_kernel<<<block_count(ctx->n_real), kThreads>>>(
      psi, ctx->n_real, 1.0 / static_cast<double>(ctx->n_global));
  if (!kernel_ok("distributed inverse FFT scale kernel")) return 1;
  if (measure_timing != 0 && !synchronize_for_timing(ctx, &stamps[5])) return 1;

  local_half_step_kernel<<<block_count(ctx->n_real), kThreads>>>(
      psi, ctx->potential, ctx->n_real, local_tau, g, imaginary_time);
  if (!kernel_ok("second local half-step kernel")) return 1;
  if (measure_timing != 0 && !synchronize_for_timing(ctx, &stamps[6])) return 1;

  if (imaginary_time != 0 &&
      !normalize_wavefunction(ctx, target_norm, volume_element)) {
    return 1;
  }
  if (measure_timing != 0) {
    if (!synchronize_for_timing(ctx, &stamps[7])) return 1;
    *nonlinear_seconds = (stamps[1] - stamps[0]) +
                         (stamps[6] - stamps[5]);
    *fft_seconds = (stamps[2] - stamps[1]) + (stamps[4] - stamps[3]);
    *kinetic_seconds = stamps[3] - stamps[2];
    *other_seconds = (stamps[5] - stamps[4]) + (stamps[7] - stamps[6]);
    *total_seconds = stamps[7] - stamps[0];
  }
  return 0;
}

GP3D_CUFFTMP_API int gp3d_cufftmp_argle_step(
    void* handle, double alpha, double reaction_constant,
    double nonlinear_coefficient, double potential_scale, double dtau,
    double velocity_amplitude, int compute_metrics, double* max_delta_rate,
    double* mean_density, double* min_density) {
  auto* ctx = static_cast<GpuContext*>(handle);
  if (ctx == nullptr || !ctx->has_argle || alpha <= 0.0 || dtau <= 0.0) {
    set_error("invalid distributed CUDA ARGLE context or parameters");
    return 1;
  }
  if (!cuda_ok(cudaSetDevice(ctx->device), "select CUDA device for ARGLE")) {
    return 1;
  }

  auto* psi = descriptor_data(ctx->psi_desc);
  auto* spectral = descriptor_data(ctx->work_desc);
  auto* old_psi = descriptor_data(ctx->old_psi_desc);
  auto* grad_x = descriptor_data(ctx->grad_x_desc);
  auto* grad_y = descriptor_data(ctx->grad_y_desc);
  auto* rhs = descriptor_data(ctx->rhs_desc);
  const size_t real_bytes = ctx->n_real * sizeof(cufftDoubleComplex);

  if (!cuda_ok(cudaMemcpy(old_psi, psi, real_bytes, cudaMemcpyDeviceToDevice),
               "copy distributed ARGLE old psi") ||
      !cuda_ok(cudaMemcpy(spectral, old_psi, real_bytes,
                          cudaMemcpyDeviceToDevice),
               "copy ARGLE psi into FFT work buffer")) {
    return 1;
  }
  if (!forward_fft(ctx, spectral, "ARGLE forward distributed transform")) {
    return 1;
  }

  derivative_shuffled_kernel<<<block_count(ctx->n_spectral), kThreads>>>(
      spectral, grad_x, ctx->kx, ctx->ky, ctx->kz, ctx->nx,
      ctx->local_ny, ctx->y_start, ctx->n_spectral, 0);
  if (!kernel_ok("ARGLE shuffled x-derivative kernel") ||
      !inverse_fft(ctx, grad_x, "ARGLE x-derivative inverse transform")) {
    return 1;
  }
  scale_kernel<<<block_count(ctx->n_real), kThreads>>>(
      grad_x, ctx->n_real, 1.0 / static_cast<double>(ctx->n_global));
  if (!kernel_ok("ARGLE x-derivative scale kernel")) return 1;

  derivative_shuffled_kernel<<<block_count(ctx->n_spectral), kThreads>>>(
      spectral, grad_y, ctx->kx, ctx->ky, ctx->kz, ctx->nx,
      ctx->local_ny, ctx->y_start, ctx->n_spectral, 1);
  if (!kernel_ok("ARGLE shuffled y-derivative kernel") ||
      !inverse_fft(ctx, grad_y, "ARGLE y-derivative inverse transform")) {
    return 1;
  }
  scale_kernel<<<block_count(ctx->n_real), kThreads>>>(
      grad_y, ctx->n_real, 1.0 / static_cast<double>(ctx->n_global));
  if (!kernel_ok("ARGLE y-derivative scale kernel")) return 1;

  argle_rhs_natural_kernel<<<block_count(ctx->n_real), kThreads>>>(
      old_psi, grad_x, grad_y, rhs, ctx->potential, ctx->nx, ctx->ny,
      ctx->local_nz, ctx->z_start, ctx->nz, ctx->n_real,
      reaction_constant, nonlinear_coefficient, potential_scale, alpha,
      velocity_amplitude);
  if (!kernel_ok("distributed ARGLE rhs kernel") ||
      !forward_fft(ctx, rhs, "ARGLE rhs forward distributed transform")) {
    return 1;
  }

  argle_update_shuffled_kernel<<<block_count(ctx->n_spectral), kThreads>>>(
      spectral, rhs, ctx->kx, ctx->ky, ctx->kz, ctx->nx, ctx->local_ny,
      ctx->y_start, ctx->n_spectral, alpha, dtau);
  if (!kernel_ok("distributed ARGLE semi-implicit update kernel") ||
      !inverse_fft(ctx, spectral, "ARGLE final inverse distributed transform")) {
    return 1;
  }
  scale_kernel<<<block_count(ctx->n_real), kThreads>>>(
      spectral, ctx->n_real, 1.0 / static_cast<double>(ctx->n_global));
  if (!kernel_ok("ARGLE final distributed scale kernel")) return 1;

  // The work descriptor now contains the new natural-layout wave function.
  std::swap(ctx->psi_desc, ctx->work_desc);
  psi = descriptor_data(ctx->psi_desc);

  *max_delta_rate = 0.0;
  *mean_density = 0.0;
  *min_density = 0.0;
  if (compute_metrics != 0) {
    const double initial[3] = {0.0, DBL_MAX, 0.0};
    if (!cuda_ok(cudaMemcpy(ctx->reductions, initial, sizeof(initial),
                            cudaMemcpyHostToDevice),
                 "initialize local ARGLE metrics")) {
      return 1;
    }
    argle_metrics_kernel<<<block_count(ctx->n_real), kThreads>>>(
        psi, old_psi, ctx->n_real, ctx->reductions);
    if (!kernel_ok("local ARGLE metrics kernel")) return 1;
    double local[3]{};
    if (!cuda_ok(cudaMemcpy(local, ctx->reductions, sizeof(local),
                            cudaMemcpyDeviceToHost),
                 "copy local ARGLE metrics")) {
      return 1;
    }
    double global_sum = 0.0;
    double global_min = 0.0;
    double global_max = 0.0;
    if (!allreduce_values(ctx, &local[0], &global_sum, 1, MPI_SUM,
                          "allreduce ARGLE density sum") ||
        !allreduce_values(ctx, &local[1], &global_min, 1, MPI_MIN,
                          "allreduce ARGLE minimum density") ||
        !allreduce_values(ctx, &local[2], &global_max, 1, MPI_MAX,
                          "allreduce ARGLE maximum change")) {
      return 1;
    }
    *mean_density = global_sum / static_cast<double>(ctx->n_global);
    *min_density = global_min;
    *max_delta_rate = global_max / dtau;
  }
  return 0;
}

GP3D_CUFFTMP_API int gp3d_cufftmp_diagnostics(
    void* handle, double hbar, double mass, double g, double volume_element,
    double* norm, double* energy) {
  auto* ctx = static_cast<GpuContext*>(handle);
  if (ctx == nullptr || mass <= 0.0) {
    set_error("invalid distributed CUDA diagnostic context or parameters");
    return 1;
  }
  if (!cuda_ok(cudaSetDevice(ctx->device),
               "select CUDA device for diagnostics") ||
      !cuda_ok(cudaMemset(ctx->reductions, 0, 3 * sizeof(double)),
               "clear local diagnostic reductions")) {
    return 1;
  }

  auto* psi = descriptor_data(ctx->psi_desc);
  auto* spectral = descriptor_data(ctx->work_desc);
  local_diagnostics_kernel<<<block_count(ctx->n_real), kThreads>>>(
      psi, ctx->potential, ctx->n_real, g, ctx->reductions);
  if (!kernel_ok("local diagnostic kernel") ||
      !cuda_ok(cudaMemcpy(spectral, psi,
                          ctx->n_real * sizeof(cufftDoubleComplex),
                          cudaMemcpyDeviceToDevice),
               "copy psi for distributed spectral diagnostics") ||
      !forward_fft(ctx, spectral, "diagnostic forward distributed transform")) {
    return 1;
  }
  spectral_diagnostics_shuffled_kernel<<<block_count(ctx->n_spectral),
                                          kThreads>>>(
      spectral, ctx->kx, ctx->ky, ctx->kz, ctx->nx, ctx->local_ny,
      ctx->y_start, ctx->n_spectral, ctx->reductions);
  if (!kernel_ok("shuffled spectral diagnostic kernel")) return 1;

  double local[3]{};
  if (!cuda_ok(cudaMemcpy(local, ctx->reductions, sizeof(local),
                          cudaMemcpyDeviceToHost),
               "copy local diagnostic reductions")) {
    return 1;
  }
  double global[3]{};
  if (!allreduce_values(ctx, local, global, 3, MPI_SUM,
                        "allreduce GPE diagnostics")) {
    return 1;
  }
  *norm = global[0] * volume_element;
  const double spectral_scale =
      volume_element / static_cast<double>(ctx->n_global);
  *energy = global[1] * volume_element +
            0.5 * hbar * hbar / mass * global[2] * spectral_scale;
  return 0;
}

GP3D_CUFFTMP_API void gp3d_cufftmp_destroy(void* handle) {
  destroy_context(static_cast<GpuContext*>(handle));
}

GP3D_CUFFTMP_API void gp3d_cufftmp_get_last_error(char* buffer,
                                                   int buffer_size) {
  if (buffer == nullptr || buffer_size <= 0) return;
  const size_t count = static_cast<size_t>(buffer_size - 1);
  std::strncpy(buffer, g_last_error, count);
  buffer[count] = '\0';
}

GP3D_CUFFTMP_API void gp3d_cufftmp_abort(int fortran_comm, int error_code) {
  const MPI_Comm comm =
      MPI_Comm_f2c(static_cast<MPI_Fint>(fortran_comm));
  MPI_Abort(comm, error_code == 0 ? 1 : error_code);
}
