// Fortran側から呼び出すCUDA/cuFFT実装。
// psi、potential、波数配列をGPUへ常駐させ、局所項・FFT・運動項・ARGLE・診断量を実行する。
#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>

#if defined(_WIN32)
#define GP3D_CUDA_API extern "C" __declspec(dllexport)
#else
#define GP3D_CUDA_API extern "C"
#endif

namespace {

constexpr int kThreads = 256;
constexpr double kPi = 3.141592653589793238462643383279502884;
char g_last_error[512] = "";

// 1ケースで使う全GPUメモリとcuFFT planを所有し、時間ループ中の再確保を避ける。
struct GpuContext {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  int device = 0;
  size_t n = 0;
  bool has_argle = false;
  bool plan_created = false;
  int events_created = 0;
  cufftHandle plan{};
  cufftDoubleComplex* psi = nullptr;
  cufftDoubleComplex* spectral = nullptr;
  cufftDoubleComplex* old_psi = nullptr;
  cufftDoubleComplex* grad_x = nullptr;
  cufftDoubleComplex* grad_y = nullptr;
  cufftDoubleComplex* rhs = nullptr;
  double* potential = nullptr;
  double* kx = nullptr;
  double* ky = nullptr;
  double* kz = nullptr;
  double* reductions = nullptr;
  cudaEvent_t events[8]{};
};

void set_error(const char* operation, cudaError_t error) {
  std::snprintf(g_last_error, sizeof(g_last_error), "%s: %s", operation,
                cudaGetErrorString(error));
}

void set_error(const char* operation, cufftResult error) {
  std::snprintf(g_last_error, sizeof(g_last_error), "%s: cuFFT error %d", operation,
                static_cast<int>(error));
}

void set_error(const char* message) {
  std::snprintf(g_last_error, sizeof(g_last_error), "%s", message);
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

bool kernel_ok(const char* operation) {
  return cuda_ok(cudaGetLastError(), operation);
}

int block_count(size_t n) {
  return static_cast<int>((n + static_cast<size_t>(kThreads) - 1) /
                          static_cast<size_t>(kThreads));
}

void destroy_context(GpuContext* ctx) {
  if (ctx == nullptr) return;
  cudaSetDevice(ctx->device);
  for (int i = 0; i < ctx->events_created; ++i) {
    cudaEventDestroy(ctx->events[i]);
  }
  if (ctx->plan_created) cufftDestroy(ctx->plan);
  cudaFree(ctx->psi);
  cudaFree(ctx->spectral);
  cudaFree(ctx->old_psi);
  cudaFree(ctx->grad_x);
  cudaFree(ctx->grad_y);
  cudaFree(ctx->rhs);
  cudaFree(ctx->potential);
  cudaFree(ctx->kx);
  cudaFree(ctx->ky);
  cudaFree(ctx->kz);
  cudaFree(ctx->reductions);
  ctx->~GpuContext();
  std::free(ctx);
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

// 実空間でポテンシャル項と非線形項の半ステップを各格子点独立に適用する。
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

// Fourier空間で対角化された運動項を各波数モードへ適用する。
__global__ void kinetic_step_kernel(cufftDoubleComplex* spectral,
                                    const double* kx, const double* ky,
                                    const double* kz, int nx, int ny, int nz,
                                    double tau, double kinetic_coefficient,
                                    int imaginary_time) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  const size_t n = static_cast<size_t>(nx) * ny * nz;
  if (index >= n) return;
  const int i = static_cast<int>(index % nx);
  const size_t plane_index = index / nx;
  const int j = static_cast<int>(plane_index % ny);
  const int k = static_cast<int>(plane_index / ny);
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

__global__ void derivative_kernel(const cufftDoubleComplex* spectral,
                                  cufftDoubleComplex* derivative,
                                  const double* wave_numbers, int nx, int ny,
                                  int nz, int axis) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  const size_t n = static_cast<size_t>(nx) * ny * nz;
  if (index >= n) return;
  const int i = static_cast<int>(index % nx);
  const size_t plane_index = index / nx;
  const int j = static_cast<int>(plane_index % ny);
  const int k = static_cast<int>(plane_index / ny);
  const int component = axis == 0 ? i : (axis == 1 ? j : k);
  const double wave_number = wave_numbers[component];
  const cufftDoubleComplex value = spectral[index];
  derivative[index] = make_cuDoubleComplex(-wave_number * value.y,
                                           wave_number * value.x);
}

// Taylor-Green移流項を含むARGLE右辺を実空間で組み立てる。
__global__ void argle_rhs_kernel(
    const cufftDoubleComplex* old_psi, const cufftDoubleComplex* grad_x,
    const cufftDoubleComplex* grad_y, cufftDoubleComplex* rhs,
    const double* potential, int nx, int ny, int nz,
    double reaction_constant, double nonlinear_coefficient,
    double potential_scale, double alpha, double velocity_amplitude) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  const size_t n = static_cast<size_t>(nx) * ny * nz;
  if (index >= n) return;
  const int i = static_cast<int>(index % nx);
  const size_t plane_index = index / nx;
  const int j = static_cast<int>(plane_index % ny);
  const int k = static_cast<int>(plane_index / ny);

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

__global__ void argle_update_kernel(cufftDoubleComplex* spectral,
                                    const cufftDoubleComplex* rhs_spectral,
                                    const double* kx, const double* ky,
                                    const double* kz, int nx, int ny, int nz,
                                    double alpha, double dtau) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  const size_t n = static_cast<size_t>(nx) * ny * nz;
  if (index >= n) return;
  const int i = static_cast<int>(index % nx);
  const size_t plane_index = index / nx;
  const int j = static_cast<int>(plane_index % ny);
  const int k = static_cast<int>(plane_index / ny);
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

__global__ void spectral_diagnostics_kernel(
    const cufftDoubleComplex* spectral, const double* kx, const double* ky,
    const double* kz, int nx, int ny, int nz, double* reductions) {
  const size_t index = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x;
  const size_t n = static_cast<size_t>(nx) * ny * nz;
  if (index >= n) return;
  const int i = static_cast<int>(index % nx);
  const size_t plane_index = index / nx;
  const int j = static_cast<int>(plane_index % ny);
  const int k = static_cast<int>(plane_index / ny);
  const double k2 = kx[i] * kx[i] + ky[j] * ky[j] + kz[k] * kz[k];
  atomic_add_double(&reductions[2], k2 * complex_abs2(spectral[index]));
}

bool normalize_wavefunction(GpuContext* ctx, double target_norm,
                            double volume_element) {
  if (!cuda_ok(cudaMemset(ctx->reductions, 0, sizeof(double)),
               "clear norm reduction")) return false;
  norm_kernel<<<block_count(ctx->n), kThreads>>>(ctx->psi, ctx->n,
                                                 ctx->reductions);
  if (!kernel_ok("norm kernel")) return false;
  double density_sum = 0.0;
  if (!cuda_ok(cudaMemcpy(&density_sum, ctx->reductions, sizeof(double),
                          cudaMemcpyDeviceToHost),
               "copy norm reduction")) return false;
  const double current_norm = density_sum * volume_element;
  if (!(current_norm > 0.0)) {
    set_error("cannot normalize a zero GPU wave function");
    return false;
  }
  const double scale = sqrt(target_norm / current_norm);
  scale_kernel<<<block_count(ctx->n), kThreads>>>(ctx->psi, ctx->n, scale);
  return kernel_ok("normalization scale kernel");
}

bool record_event(GpuContext* ctx, int index) {
  return cuda_ok(cudaEventRecord(ctx->events[index]), "record CUDA timing event");
}

double elapsed_seconds(GpuContext* ctx, int first, int second) {
  float milliseconds = 0.0F;
  cudaEventElapsedTime(&milliseconds, ctx->events[first], ctx->events[second]);
  return static_cast<double>(milliseconds) * 1.0e-3;
}

}  // namespace

// 以下はISO_C_BINDINGから呼ばれるC ABI。例外を越境させず、失敗はstatusと文字列で返す。
GP3D_CUDA_API int gp3d_cuda_create(void** handle, int nx, int ny, int nz,
                                   int need_argle, int device) {
  g_last_error[0] = '\0';
  if (handle == nullptr || nx <= 0 || ny <= 0 || nz <= 0) {
    set_error("invalid CUDA context dimensions");
    return 1;
  }
  void* storage = std::malloc(sizeof(GpuContext));
  if (storage == nullptr) {
    set_error("failed to allocate CUDA context");
    return 1;
  }
  auto* ctx = ::new (storage) GpuContext();
  ctx->nx = nx;
  ctx->ny = ny;
  ctx->nz = nz;
  ctx->device = device;
  ctx->n = static_cast<size_t>(nx) * ny * nz;
  ctx->has_argle = need_argle != 0;
  const size_t complex_bytes = ctx->n * sizeof(cufftDoubleComplex);
  const size_t real_bytes = ctx->n * sizeof(double);

  if (!cuda_ok(cudaSetDevice(device), "select CUDA device")) {
    destroy_context(ctx);
    return 1;
  }
#define GP3D_CUDA_ALLOC(member, bytes, label)                                  \
  do {                                                                          \
    if (!cuda_ok(cudaMalloc(reinterpret_cast<void**>(&(member)), (bytes)),       \
                 (label))) {                                                    \
      destroy_context(ctx);                                                     \
      return 1;                                                                 \
    }                                                                           \
  } while (false)

  GP3D_CUDA_ALLOC(ctx->psi, complex_bytes, "allocate GPU psi");
  GP3D_CUDA_ALLOC(ctx->spectral, complex_bytes, "allocate GPU spectral buffer");
  GP3D_CUDA_ALLOC(ctx->potential, real_bytes, "allocate GPU potential");
  GP3D_CUDA_ALLOC(ctx->kx, static_cast<size_t>(nx) * sizeof(double),
                  "allocate GPU kx");
  GP3D_CUDA_ALLOC(ctx->ky, static_cast<size_t>(ny) * sizeof(double),
                  "allocate GPU ky");
  GP3D_CUDA_ALLOC(ctx->kz, static_cast<size_t>(nz) * sizeof(double),
                  "allocate GPU kz");
  GP3D_CUDA_ALLOC(ctx->reductions, 3 * sizeof(double),
                  "allocate GPU reductions");
  if (ctx->has_argle) {
    GP3D_CUDA_ALLOC(ctx->old_psi, complex_bytes, "allocate GPU ARGLE old psi");
    GP3D_CUDA_ALLOC(ctx->grad_x, complex_bytes, "allocate GPU ARGLE grad x");
    GP3D_CUDA_ALLOC(ctx->grad_y, complex_bytes, "allocate GPU ARGLE grad y");
    GP3D_CUDA_ALLOC(ctx->rhs, complex_bytes, "allocate GPU ARGLE rhs");
  }
#undef GP3D_CUDA_ALLOC

  if (!cufft_ok(cufftPlan3d(&ctx->plan, nz, ny, nx, CUFFT_Z2Z),
                "create 3D cuFFT plan")) {
    destroy_context(ctx);
    return 1;
  }
  ctx->plan_created = true;
  for (int i = 0; i < 8; ++i) {
    if (!cuda_ok(cudaEventCreate(&ctx->events[i]), "create CUDA timing event")) {
      destroy_context(ctx);
      return 1;
    }
    ++ctx->events_created;
  }
  *handle = ctx;
  return 0;
}

GP3D_CUDA_API int gp3d_cuda_upload(void* handle,
                                   const cufftDoubleComplex* psi,
                                   const double* potential, const double* kx,
                                   const double* ky, const double* kz) {
  auto* ctx = static_cast<GpuContext*>(handle);
  if (ctx == nullptr) {
    set_error("CUDA context is null during upload");
    return 1;
  }
  if (!cuda_ok(cudaSetDevice(ctx->device), "select CUDA device")) return 1;
  if (!cuda_ok(cudaMemcpy(ctx->psi, psi,
                          ctx->n * sizeof(cufftDoubleComplex),
                          cudaMemcpyHostToDevice), "upload psi")) return 1;
  if (!cuda_ok(cudaMemcpy(ctx->potential, potential,
                          ctx->n * sizeof(double), cudaMemcpyHostToDevice),
               "upload potential")) return 1;
  if (!cuda_ok(cudaMemcpy(ctx->kx, kx, ctx->nx * sizeof(double),
                          cudaMemcpyHostToDevice), "upload kx")) return 1;
  if (!cuda_ok(cudaMemcpy(ctx->ky, ky, ctx->ny * sizeof(double),
                          cudaMemcpyHostToDevice), "upload ky")) return 1;
  if (!cuda_ok(cudaMemcpy(ctx->kz, kz, ctx->nz * sizeof(double),
                          cudaMemcpyHostToDevice), "upload kz")) return 1;
  return 0;
}

GP3D_CUDA_API int gp3d_cuda_download(void* handle,
                                     cufftDoubleComplex* psi) {
  auto* ctx = static_cast<GpuContext*>(handle);
  if (ctx == nullptr) {
    set_error("CUDA context is null during download");
    return 1;
  }
  if (!cuda_ok(cudaSetDevice(ctx->device), "select CUDA device")) return 1;
  return cuda_ok(cudaMemcpy(psi, ctx->psi,
                            ctx->n * sizeof(cufftDoubleComplex),
                            cudaMemcpyDeviceToHost), "download psi") ? 0 : 1;
}

// 局所半ステップ、前進cuFFT、運動項、逆cuFFT、局所半ステップを一括実行する。
GP3D_CUDA_API int gp3d_cuda_step(
    void* handle, double dt, double hbar, double mass, double g,
    double target_norm, double volume_element, int imaginary_time,
    int measure_timing, double* nonlinear_seconds, double* fft_seconds,
    double* kinetic_seconds, double* other_seconds, double* total_seconds) {
  auto* ctx = static_cast<GpuContext*>(handle);
  if (ctx == nullptr || hbar <= 0.0 || mass <= 0.0) {
    set_error("invalid CUDA split-step context or parameters");
    return 1;
  }
  if (!cuda_ok(cudaSetDevice(ctx->device), "select CUDA device")) return 1;
  *nonlinear_seconds = 0.0;
  *fft_seconds = 0.0;
  *kinetic_seconds = 0.0;
  *other_seconds = 0.0;
  *total_seconds = 0.0;
  if (measure_timing != 0 && !record_event(ctx, 0)) return 1;

  const double local_tau = 0.5 * dt / hbar;
  local_half_step_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->psi, ctx->potential, ctx->n, local_tau, g, imaginary_time);
  if (!kernel_ok("first local half-step kernel")) return 1;
  if (measure_timing != 0 && !record_event(ctx, 1)) return 1;

  if (!cufft_ok(cufftExecZ2Z(ctx->plan, ctx->psi, ctx->spectral,
                             CUFFT_FORWARD), "forward cuFFT")) return 1;
  if (measure_timing != 0 && !record_event(ctx, 2)) return 1;

  kinetic_step_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->spectral, ctx->kx, ctx->ky, ctx->kz, ctx->nx, ctx->ny, ctx->nz,
      dt / hbar, 0.5 * hbar * hbar / mass, imaginary_time);
  if (!kernel_ok("kinetic spectral kernel")) return 1;
  if (measure_timing != 0 && !record_event(ctx, 3)) return 1;

  if (!cufft_ok(cufftExecZ2Z(ctx->plan, ctx->spectral, ctx->psi,
                             CUFFT_INVERSE), "inverse cuFFT")) return 1;
  if (measure_timing != 0 && !record_event(ctx, 4)) return 1;

  scale_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->psi, ctx->n, 1.0 / static_cast<double>(ctx->n));
  if (!kernel_ok("inverse FFT scale kernel")) return 1;
  if (measure_timing != 0 && !record_event(ctx, 5)) return 1;

  local_half_step_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->psi, ctx->potential, ctx->n, local_tau, g, imaginary_time);
  if (!kernel_ok("second local half-step kernel")) return 1;
  if (measure_timing != 0 && !record_event(ctx, 6)) return 1;

  if (imaginary_time != 0 &&
      !normalize_wavefunction(ctx, target_norm, volume_element)) return 1;
  if (measure_timing != 0) {
    if (!record_event(ctx, 7)) return 1;
    if (!cuda_ok(cudaEventSynchronize(ctx->events[7]),
                 "synchronize CUDA timing event")) return 1;
    *nonlinear_seconds = elapsed_seconds(ctx, 0, 1) +
                         elapsed_seconds(ctx, 5, 6);
    *fft_seconds = elapsed_seconds(ctx, 1, 2) + elapsed_seconds(ctx, 3, 4);
    *kinetic_seconds = elapsed_seconds(ctx, 2, 3);
    *other_seconds = elapsed_seconds(ctx, 4, 5) + elapsed_seconds(ctx, 6, 7);
    *total_seconds = elapsed_seconds(ctx, 0, 7);
  }
  return 0;
}

// 半陰的ARGLEを1ステップ進め、要求時だけ収束指標をGPU上で集約する。
GP3D_CUDA_API int gp3d_cuda_argle_step(
    void* handle, double alpha, double reaction_constant,
    double nonlinear_coefficient, double potential_scale, double dtau,
    double velocity_amplitude, int compute_metrics, double* max_delta_rate,
    double* mean_density, double* min_density) {
  auto* ctx = static_cast<GpuContext*>(handle);
  if (ctx == nullptr || !ctx->has_argle || alpha <= 0.0 || dtau <= 0.0) {
    set_error("invalid CUDA ARGLE context or parameters");
    return 1;
  }
  if (!cuda_ok(cudaSetDevice(ctx->device), "select CUDA device")) return 1;
  const size_t complex_bytes = ctx->n * sizeof(cufftDoubleComplex);
  if (!cuda_ok(cudaMemcpy(ctx->old_psi, ctx->psi, complex_bytes,
                          cudaMemcpyDeviceToDevice), "copy ARGLE old psi")) return 1;
  if (!cufft_ok(cufftExecZ2Z(ctx->plan, ctx->old_psi, ctx->spectral,
                             CUFFT_FORWARD), "ARGLE forward cuFFT")) return 1;

  derivative_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->spectral, ctx->grad_x, ctx->kx, ctx->nx, ctx->ny, ctx->nz, 0);
  if (!kernel_ok("ARGLE x derivative kernel")) return 1;
  if (!cufft_ok(cufftExecZ2Z(ctx->plan, ctx->grad_x, ctx->grad_x,
                             CUFFT_INVERSE), "ARGLE x inverse cuFFT")) return 1;
  scale_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->grad_x, ctx->n, 1.0 / static_cast<double>(ctx->n));
  if (!kernel_ok("ARGLE x derivative scale kernel")) return 1;

  derivative_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->spectral, ctx->grad_y, ctx->ky, ctx->nx, ctx->ny, ctx->nz, 1);
  if (!kernel_ok("ARGLE y derivative kernel")) return 1;
  if (!cufft_ok(cufftExecZ2Z(ctx->plan, ctx->grad_y, ctx->grad_y,
                             CUFFT_INVERSE), "ARGLE y inverse cuFFT")) return 1;
  scale_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->grad_y, ctx->n, 1.0 / static_cast<double>(ctx->n));
  if (!kernel_ok("ARGLE y derivative scale kernel")) return 1;

  argle_rhs_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->old_psi, ctx->grad_x, ctx->grad_y, ctx->rhs, ctx->potential,
      ctx->nx, ctx->ny, ctx->nz, reaction_constant, nonlinear_coefficient,
      potential_scale, alpha, velocity_amplitude);
  if (!kernel_ok("ARGLE rhs kernel")) return 1;
  if (!cufft_ok(cufftExecZ2Z(ctx->plan, ctx->rhs, ctx->rhs,
                             CUFFT_FORWARD), "ARGLE rhs forward cuFFT")) return 1;

  argle_update_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->spectral, ctx->rhs, ctx->kx, ctx->ky, ctx->kz, ctx->nx, ctx->ny,
      ctx->nz, alpha, dtau);
  if (!kernel_ok("ARGLE semi-implicit update kernel")) return 1;
  if (!cufft_ok(cufftExecZ2Z(ctx->plan, ctx->spectral, ctx->psi,
                             CUFFT_INVERSE), "ARGLE final inverse cuFFT")) return 1;
  scale_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->psi, ctx->n, 1.0 / static_cast<double>(ctx->n));
  if (!kernel_ok("ARGLE final scale kernel")) return 1;

  *max_delta_rate = 0.0;
  *mean_density = 0.0;
  *min_density = 0.0;
  if (compute_metrics != 0) {
    const double initial[3] = {0.0, DBL_MAX, 0.0};
    if (!cuda_ok(cudaMemcpy(ctx->reductions, initial, sizeof(initial),
                            cudaMemcpyHostToDevice),
                 "initialize ARGLE metrics")) return 1;
    argle_metrics_kernel<<<block_count(ctx->n), kThreads>>>(
        ctx->psi, ctx->old_psi, ctx->n, ctx->reductions);
    if (!kernel_ok("ARGLE metrics kernel")) return 1;
    double metrics[3];
    if (!cuda_ok(cudaMemcpy(metrics, ctx->reductions, sizeof(metrics),
                            cudaMemcpyDeviceToHost),
                 "copy ARGLE metrics")) return 1;
    *mean_density = metrics[0] / static_cast<double>(ctx->n);
    *min_density = metrics[1];
    *max_delta_rate = metrics[2] / dtau;
  }
  return 0;
}

GP3D_CUDA_API int gp3d_cuda_diagnostics(void* handle, double hbar,
                                        double mass, double g,
                                        double volume_element, double* norm,
                                        double* energy) {
  auto* ctx = static_cast<GpuContext*>(handle);
  if (ctx == nullptr || mass <= 0.0) {
    set_error("invalid CUDA diagnostic context or parameters");
    return 1;
  }
  if (!cuda_ok(cudaSetDevice(ctx->device), "select CUDA device")) return 1;
  if (!cuda_ok(cudaMemset(ctx->reductions, 0, 3 * sizeof(double)),
               "clear diagnostic reductions")) return 1;
  local_diagnostics_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->psi, ctx->potential, ctx->n, g, ctx->reductions);
  if (!kernel_ok("local diagnostic kernel")) return 1;
  if (!cufft_ok(cufftExecZ2Z(ctx->plan, ctx->psi, ctx->spectral,
                             CUFFT_FORWARD), "diagnostic forward cuFFT")) return 1;
  spectral_diagnostics_kernel<<<block_count(ctx->n), kThreads>>>(
      ctx->spectral, ctx->kx, ctx->ky, ctx->kz, ctx->nx, ctx->ny, ctx->nz,
      ctx->reductions);
  if (!kernel_ok("spectral diagnostic kernel")) return 1;
  double values[3];
  if (!cuda_ok(cudaMemcpy(values, ctx->reductions, sizeof(values),
                          cudaMemcpyDeviceToHost),
               "copy diagnostic reductions")) return 1;
  *norm = values[0] * volume_element;
  const double spectral_scale = volume_element / static_cast<double>(ctx->n);
  *energy = values[1] * volume_element +
            0.5 * hbar * hbar / mass * values[2] * spectral_scale;
  return 0;
}

GP3D_CUDA_API void gp3d_cuda_destroy(void* handle) {
  destroy_context(static_cast<GpuContext*>(handle));
}

GP3D_CUDA_API void gp3d_cuda_get_last_error(char* buffer, int buffer_size) {
  if (buffer == nullptr || buffer_size <= 0) return;
  const size_t count = static_cast<size_t>(buffer_size - 1);
  std::strncpy(buffer, g_last_error, count);
  buffer[count] = '\0';
}
