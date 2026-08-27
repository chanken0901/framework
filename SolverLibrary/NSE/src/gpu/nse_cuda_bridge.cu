#include <cuda_runtime.h>
#include <cub/device/device_reduce.cuh>
#if defined(NSE_FORCING_CUFFTMP) || defined(NSE_INIT_CUFFTMP)
#include <cufftMp.h>
#include <mpi.h>
#elif defined(NSE_FORCING_CUFFT) || defined(NSE_INIT_CUFFT)
#include <cufft.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <new>
#include <string>

#if defined(_WIN32)
#define NSE_CUDA_EXPORT extern "C" __declspec(dllexport)
#else
#define NSE_CUDA_EXPORT extern "C"
#endif

namespace {

constexpr int convective_keep2 = 1;
constexpr int convective_keep6 = 2;
constexpr int convective_weno5z_roe = 3;
constexpr int convective_hybrid = 4;
constexpr int sensor_ducros_pressure = 1;
constexpr int halo_direction_y = 1;
constexpr int halo_direction_z = 2;
constexpr int halo_side_low = -1;
constexpr int halo_side_high = 1;

thread_local std::string last_error;

struct GridView {
  int nx;
  int ny;
  int nz;
  int nghost;
  int nx_total;
  int ny_total;
  int nz_total;
  std::size_t cell_count;
};

struct NseCudaContext {
  GridView grid{};
  int nvar = 0;
  int device = 0;
  bool distributed_y = false;
  bool distributed_z = false;
  int global_ny = 0;
  int global_nz = 0;
  int global_y_start = 0;
  int global_z_start = 0;
  int convective_scheme = convective_keep6;
  int hybrid_smooth_scheme = convective_keep6;
  int hybrid_shock_scheme = convective_weno5z_roe;
  int hybrid_sensor = sensor_ducros_pressure;
  double hybrid_sensor_onset = 0.01;
  double hybrid_sensor_full = 0.10;
  double gamma = 0.0;
  double cfl = 0.0;
  double small_rho = 0.0;
  double small_p = 0.0;
  double reynolds = 0.0;
  double prandtl = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  bool viscous_enabled = false;
  bool forcing_enabled = false;
  int forcing_spectrum = 0;
  int forcing_report_interval = 0;
  int forcing_evaluations = 0;
  double forcing_k_cutoff = 0.0;
  double forcing_target_dissipation = 0.0;
  double forcing_dilatational_ratio = 0.0;
  double forcing_denominator_floor = 0.0;
  double forcing_max_coefficient = 0.0;
  std::size_t physical_count = 0;
  std::size_t state_bytes = 0;
  double* q = nullptr;
  double* q0 = nullptr;
  double* rhs = nullptr;
  double* primitive = nullptr;
  double* speed = nullptr;
  double* max_speed = nullptr;
  double* halo_buffer = nullptr;
  std::size_t halo_buffer_count = 0;
  void* reduce_storage = nullptr;
  std::size_t reduce_storage_bytes = 0;
#if defined(NSE_FORCING_CUFFT)
  cufftHandle forcing_plan = 0;
  cufftDoubleComplex* forcing_spectral = nullptr;
  cufftDoubleComplex* forcing_phi = nullptr;
  double* forcing_energy_d = nullptr;
  double* forcing_sum = nullptr;
#endif
#if defined(NSE_FORCING_CUFFTMP)
  MPI_Comm forcing_comm = MPI_COMM_NULL;
  int forcing_rank = 0;
  cufftHandle forcing_plan = 0;
  cudaLibXtDesc* forcing_physical[3] = {nullptr, nullptr, nullptr};
  cudaLibXtDesc* forcing_spectral[3] = {nullptr, nullptr, nullptr};
  cufftDoubleComplex* forcing_phi = nullptr;
  double* forcing_energy_d = nullptr;
  double* forcing_sum = nullptr;
  bool forcing_distributed_ready = false;
#endif
};

#if defined(NSE_INIT_CUFFTMP)
struct HitCufftMpContext {
  MPI_Comm comm = MPI_COMM_NULL;
  cufftHandle plan = 0;
  cudaLibXtDesc* physical = nullptr;
  cudaLibXtDesc* spectral = nullptr;
  std::size_t local_count = 0;
};

HitCufftMpContext hit_cufftmp;
#endif

void set_error(const std::string& operation, cudaError_t status) {
  last_error = operation + ": " + cudaGetErrorString(status);
}

void set_error(const std::string& message) {
  last_error = message;
}

bool check_cuda(cudaError_t status, const char* operation) {
  if (status == cudaSuccess) {
    return true;
  }
  set_error(operation, status);
  return false;
}

#if defined(NSE_FORCING_CUFFT) || defined(NSE_INIT_CUFFT) || \
    defined(NSE_FORCING_CUFFTMP) || defined(NSE_INIT_CUFFTMP)
bool check_cufft(cufftResult status, const char* operation) {
  if (status == CUFFT_SUCCESS) {
    return true;
  }
  last_error = std::string(operation) + ": cuFFT status "
      + std::to_string(static_cast<int>(status));
  return false;
}
#endif

#if defined(NSE_FORCING_CUFFTMP) || defined(NSE_INIT_CUFFTMP)
cufftDoubleComplex* descriptor_data(cudaLibXtDesc* descriptor) {
  if (descriptor == nullptr || descriptor->descriptor == nullptr
      || descriptor->descriptor->data[0] == nullptr) {
    set_error("cuFFTMp descriptor has no local device allocation");
    return nullptr;
  }
  return static_cast<cufftDoubleComplex*>(descriptor->descriptor->data[0]);
}
#endif

void release_context(NseCudaContext* context) {
  if (context == nullptr) {
    return;
  }
  cudaSetDevice(context->device);
#if defined(NSE_FORCING_CUFFT)
  if (context->forcing_plan != 0) {
    cufftDestroy(context->forcing_plan);
  }
  cudaFree(context->forcing_sum);
  cudaFree(context->forcing_energy_d);
  cudaFree(context->forcing_phi);
  cudaFree(context->forcing_spectral);
#endif
#if defined(NSE_FORCING_CUFFTMP)
  for (int component = 0; component < 3; ++component) {
    if (context->forcing_spectral[component] != nullptr) {
      cufftXtFree(context->forcing_spectral[component]);
    }
    if (context->forcing_physical[component] != nullptr) {
      cufftXtFree(context->forcing_physical[component]);
    }
  }
  if (context->forcing_plan != 0) cufftDestroy(context->forcing_plan);
  cudaFree(context->forcing_sum);
  cudaFree(context->forcing_energy_d);
  cudaFree(context->forcing_phi);
  if (context->forcing_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&context->forcing_comm);
  }
#endif
  cudaFree(context->reduce_storage);
  cudaFree(context->halo_buffer);
  cudaFree(context->max_speed);
  cudaFree(context->speed);
  cudaFree(context->primitive);
  cudaFree(context->rhs);
  cudaFree(context->q0);
  cudaFree(context->q);
  delete context;
}

__host__ __device__ inline std::size_t cell_index(
    const GridView& grid, int x, int y, int z) {
  return static_cast<std::size_t>(x)
      + static_cast<std::size_t>(grid.nx_total)
          * (static_cast<std::size_t>(y)
              + static_cast<std::size_t>(grid.ny_total)
                  * static_cast<std::size_t>(z));
}

__host__ __device__ inline std::size_t state_index(
    const GridView& grid, int x, int y, int z, int variable) {
  return cell_index(grid, x, y, z)
      + grid.cell_count * static_cast<std::size_t>(variable);
}

__device__ inline int positive_mod(int value, int modulus) {
  const int result = value % modulus;
  return result < 0 ? result + modulus : result;
}

__global__ void periodic_halo_kernel(
    double* q, GridView grid, int nvar, int wrap_y, int wrap_z) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= grid.cell_count) {
    return;
  }

  const int x = static_cast<int>(linear % grid.nx_total);
  const std::size_t yz = linear / grid.nx_total;
  const int y = static_cast<int>(yz % grid.ny_total);
  const int z = static_cast<int>(yz / grid.ny_total);
  const bool x_halo = x < grid.nghost || x >= grid.nghost + grid.nx;
  const bool y_halo = y < grid.nghost || y >= grid.nghost + grid.ny;
  const bool z_halo = z < grid.nghost || z >= grid.nghost + grid.nz;
  const bool interior = !x_halo && !y_halo && !z_halo;
  if (interior) {
    return;
  }
  if ((y_halo && wrap_y == 0) || (z_halo && wrap_z == 0)) {
    return;
  }

  const int wrapped_x =
      grid.nghost + positive_mod(x - grid.nghost, grid.nx);
  const int wrapped_y = wrap_y != 0
      ? grid.nghost + positive_mod(y - grid.nghost, grid.ny) : y;
  const int wrapped_z = wrap_z != 0
      ? grid.nghost + positive_mod(z - grid.nghost, grid.nz) : z;
  for (int variable = 0; variable < nvar; ++variable) {
    q[state_index(grid, x, y, z, variable)] =
        q[state_index(grid, wrapped_x, wrapped_y, wrapped_z, variable)];
  }
}

__global__ void pack_y_halo_kernel(
    const double* q, double* buffer, GridView grid, int nvar, int side,
    std::size_t count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= count) return;

  std::size_t value = linear;
  const int x = static_cast<int>(value % grid.nx_total);
  value /= grid.nx_total;
  const int layer = static_cast<int>(value % grid.nghost);
  value /= grid.nghost;
  const int z = static_cast<int>(value % grid.nz_total);
  const int variable = static_cast<int>(value / grid.nz_total);
  const int y = side == halo_side_low
      ? grid.nghost + layer
      : grid.nghost + grid.ny - grid.nghost + layer;
  buffer[linear] = q[state_index(grid, x, y, z, variable)];
}

__global__ void unpack_y_halo_kernel(
    double* q, const double* buffer, GridView grid, int nvar, int side,
    std::size_t count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= count) return;

  std::size_t value = linear;
  const int x = static_cast<int>(value % grid.nx_total);
  value /= grid.nx_total;
  const int layer = static_cast<int>(value % grid.nghost);
  value /= grid.nghost;
  const int z = static_cast<int>(value % grid.nz_total);
  const int variable = static_cast<int>(value / grid.nz_total);
  const int y = side == halo_side_low
      ? layer : grid.nghost + grid.ny + layer;
  q[state_index(grid, x, y, z, variable)] = buffer[linear];
}

__global__ void pack_z_halo_kernel(
    const double* q, double* buffer, GridView grid, int nvar, int side,
    std::size_t count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= count) return;

  std::size_t value = linear;
  const int x = static_cast<int>(value % grid.nx_total);
  value /= grid.nx_total;
  const int y = static_cast<int>(value % grid.ny_total);
  value /= grid.ny_total;
  const int layer = static_cast<int>(value % grid.nghost);
  const int variable = static_cast<int>(value / grid.nghost);
  const int z = side == halo_side_low
      ? grid.nghost + layer
      : grid.nghost + grid.nz - grid.nghost + layer;
  buffer[linear] = q[state_index(grid, x, y, z, variable)];
}

__global__ void unpack_z_halo_kernel(
    double* q, const double* buffer, GridView grid, int nvar, int side,
    std::size_t count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= count) return;

  std::size_t value = linear;
  const int x = static_cast<int>(value % grid.nx_total);
  value /= grid.nx_total;
  const int y = static_cast<int>(value % grid.ny_total);
  value /= grid.ny_total;
  const int layer = static_cast<int>(value % grid.nghost);
  const int variable = static_cast<int>(value / grid.nghost);
  const int z = side == halo_side_low
      ? layer : grid.nghost + grid.nz + layer;
  q[state_index(grid, x, y, z, variable)] = buffer[linear];
}

__device__ inline void primitive_state(
    const double* q,
    const GridView& grid,
    std::size_t cell,
    double gamma,
    double small_rho,
    double small_p,
    double& rho,
    double& u,
    double& v,
    double& w,
    double& pressure) {
  rho = fmax(q[cell], small_rho);
  u = q[cell + grid.cell_count] / rho;
  v = q[cell + 2 * grid.cell_count] / rho;
  w = q[cell + 3 * grid.cell_count] / rho;
  const double energy = q[cell + 4 * grid.cell_count];
  pressure = fmax(
      (gamma - 1.0)
          * (energy - 0.5 * rho * (u * u + v * v + w * w)),
      small_p);
}

__global__ void primitive_kernel(
    const double* q,
    double* primitive,
    GridView grid,
    double gamma,
    double small_rho,
    double small_p) {
  const std::size_t cell =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (cell >= grid.cell_count) {
    return;
  }

  double rho;
  double u;
  double v;
  double w;
  double pressure;
  primitive_state(
      q, grid, cell, gamma, small_rho, small_p,
      rho, u, v, w, pressure);
  primitive[cell] = u;
  primitive[cell + grid.cell_count] = v;
  primitive[cell + 2 * grid.cell_count] = w;
  primitive[cell + 3 * grid.cell_count] = pressure / rho;
}

__device__ inline double scalar_value(
    const double* field, const GridView& grid, int i, int j, int k) {
  return field[cell_index(grid, i, j, k)];
}

__device__ inline double first_derivative(
    const double* field,
    const GridView& grid,
    int i,
    int j,
    int k,
    int direction,
    double inverse_spacing) {
  const int di = direction == 0 ? 1 : 0;
  const int dj = direction == 1 ? 1 : 0;
  const int dk = direction == 2 ? 1 : 0;
  return inverse_spacing / 60.0 * (
      -scalar_value(field, grid, i - 3 * di, j - 3 * dj, k - 3 * dk)
      + 9.0 * scalar_value(
          field, grid, i - 2 * di, j - 2 * dj, k - 2 * dk)
      - 45.0 * scalar_value(
          field, grid, i - di, j - dj, k - dk)
      + 45.0 * scalar_value(
          field, grid, i + di, j + dj, k + dk)
      - 9.0 * scalar_value(
          field, grid, i + 2 * di, j + 2 * dj, k + 2 * dk)
      + scalar_value(
          field, grid, i + 3 * di, j + 3 * dj, k + 3 * dk));
}

__device__ inline double second_derivative(
    const double* field,
    const GridView& grid,
    int i,
    int j,
    int k,
    int direction,
    double inverse_spacing_square) {
  const int di = direction == 0 ? 1 : 0;
  const int dj = direction == 1 ? 1 : 0;
  const int dk = direction == 2 ? 1 : 0;
  return inverse_spacing_square / 180.0 * (
      2.0 * scalar_value(
          field, grid, i - 3 * di, j - 3 * dj, k - 3 * dk)
      - 27.0 * scalar_value(
          field, grid, i - 2 * di, j - 2 * dj, k - 2 * dk)
      + 270.0 * scalar_value(
          field, grid, i - di, j - dj, k - dk)
      - 490.0 * scalar_value(field, grid, i, j, k)
      + 270.0 * scalar_value(
          field, grid, i + di, j + dj, k + dk)
      - 27.0 * scalar_value(
          field, grid, i + 2 * di, j + 2 * dj, k + 2 * dk)
      + 2.0 * scalar_value(
          field, grid, i + 3 * di, j + 3 * dj, k + 3 * dk));
}

__device__ inline double mixed_derivative(
    const double* field,
    const GridView& grid,
    int i,
    int j,
    int k,
    int first_direction,
    int second_direction,
    double inverse_spacing_product) {
  constexpr double coefficient[3] = {
      3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0};
  double value = 0.0;
  for (int b = 1; b <= 3; ++b) {
    const int bi = second_direction == 0 ? b : 0;
    const int bj = second_direction == 1 ? b : 0;
    const int bk = second_direction == 2 ? b : 0;
    for (int a = 1; a <= 3; ++a) {
      const int ai = first_direction == 0 ? a : 0;
      const int aj = first_direction == 1 ? a : 0;
      const int ak = first_direction == 2 ? a : 0;
      value += coefficient[a - 1] * coefficient[b - 1] * (
          scalar_value(field, grid, i + ai + bi, j + aj + bj, k + ak + bk)
          - scalar_value(
              field, grid, i + ai - bi, j + aj - bj, k + ak - bk)
          - scalar_value(
              field, grid, i - ai + bi, j - aj + bj, k - ak + bk)
          + scalar_value(
              field, grid, i - ai - bi, j - aj - bj, k - ak - bk));
    }
  }
  return value * inverse_spacing_product;
}

__device__ inline void keep_flux(
    const double* q,
    const GridView& grid,
    std::size_t minus_cell,
    std::size_t plus_cell,
    int direction,
    double gamma,
    double small_rho,
    double small_p,
    double flux[5]) {
  double rm;
  double um;
  double vm;
  double wm;
  double pm;
  double rp;
  double up;
  double vp;
  double wp;
  double pp;
  primitive_state(
      q, grid, minus_cell, gamma, small_rho, small_p, rm, um, vm, wm, pm);
  primitive_state(
      q, grid, plus_cell, gamma, small_rho, small_p, rp, up, vp, wp, pp);

  const double rmp = 0.5 * (rm + rp);
  const double ump = 0.5 * (um + up);
  const double vmp = 0.5 * (vm + vp);
  const double wmp = 0.5 * (wm + wp);
  const double lmp =
      0.5 * (pm / rm + pp / rp) / (gamma - 1.0);
  const double normal_velocity =
      direction == 0 ? ump : (direction == 1 ? vmp : wmp);
  const double ck = rmp * normal_velocity;
  const double mxk = ck * ump;
  const double myk = ck * vmp;
  const double mzk = ck * wmp;
  const double kk = ck * 0.5 * (um * up + vm * vp + wm * wp);
  const double lk = ck * lmp;
  const double gk = 0.5 * (pm + pp);
  const double minus_normal =
      direction == 0 ? um : (direction == 1 ? vm : wm);
  const double plus_normal =
      direction == 0 ? up : (direction == 1 ? vp : wp);
  const double pk = 0.5 * (plus_normal * pm + minus_normal * pp);

  flux[0] = ck;
  flux[1] = mxk + (direction == 0 ? gk : 0.0);
  flux[2] = myk + (direction == 1 ? gk : 0.0);
  flux[3] = mzk + (direction == 2 ? gk : 0.0);
  flux[4] = kk + lk + pk;
}

__device__ inline void keep_face_flux_cuda(
    const double* q,
    const GridView& grid,
    int i,
    int j,
    int k,
    int direction,
    int keep_order,
    double gamma,
    double small_rho,
    double small_p,
    double flux[5]) {
  constexpr double central6_coefficient[3] = {
      3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0};
  const int maximum_separation = keep_order == 2 ? 1 : 3;
  for (int variable = 0; variable < 5; ++variable) {
    flux[variable] = 0.0;
  }
  for (int separation = 1; separation <= maximum_separation; ++separation) {
    const double derivative_coefficient =
        keep_order == 2 ? 0.5 : central6_coefficient[separation - 1];
    const double weight = 2.0 * derivative_coefficient;
    for (int offset = 0; offset < separation; ++offset) {
      const int im = i - (direction == 0 ? offset : 0);
      const int jm = j - (direction == 1 ? offset : 0);
      const int km = k - (direction == 2 ? offset : 0);
      const int ip = i + (direction == 0 ? separation - offset : 0);
      const int jp = j + (direction == 1 ? separation - offset : 0);
      const int kp = k + (direction == 2 ? separation - offset : 0);
      double pair_flux[5];
      keep_flux(
          q, grid, cell_index(grid, im, jm, km),
          cell_index(grid, ip, jp, kp), direction,
          gamma, small_rho, small_p, pair_flux);
      for (int variable = 0; variable < 5; ++variable) {
        flux[variable] += weight * pair_flux[variable];
      }
    }
  }
}

__global__ void wave_speed_kernel(
    const double* q,
    double* speed,
    GridView grid,
    double gamma,
    double cfl,
    double small_rho,
    double small_p,
    double reynolds,
    double prandtl,
    double inverse_dx,
    double inverse_dy,
    double inverse_dz,
    bool viscous_enabled,
    std::size_t physical_count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= physical_count) {
    return;
  }

  const int i = static_cast<int>(linear % grid.nx);
  const std::size_t jk = linear / grid.nx;
  const int j = static_cast<int>(jk % grid.ny);
  const int k = static_cast<int>(jk / grid.ny);
  const std::size_t cell = cell_index(
      grid, i + grid.nghost, j + grid.nghost, k + grid.nghost);
  double rho;
  double u;
  double v;
  double w;
  double pressure;
  primitive_state(
      q, grid, cell, gamma, small_rho, small_p,
      rho, u, v, w, pressure);
  const double sound_speed = sqrt(gamma * pressure / rho);
  const double inverse_minimum_spacing =
      fmax(inverse_dx, fmax(inverse_dy, inverse_dz));
  const double convective_rate = inverse_minimum_spacing * fmax(
      fabs(u) + sound_speed,
      fmax(fabs(v) + sound_speed, fabs(w) + sound_speed));
  double rate = convective_rate;
  if (viscous_enabled) {
    constexpr double d2_spectral_radius = 272.0 / 45.0;
    constexpr double diffusion_stability_radius = 2.0;
    const double maximum_diffusivity =
        fmax(4.0 / 3.0, gamma / prandtl) / (reynolds * rho);
    const double diffusion_rate =
        d2_spectral_radius * maximum_diffusivity
        * (inverse_dx * inverse_dx
           + inverse_dy * inverse_dy
           + inverse_dz * inverse_dz)
        / diffusion_stability_radius;
    rate = fmax(rate, cfl * diffusion_rate);
  }
  speed[linear] = rate;
}

__global__ void rhs_keep_kernel(
    const double* q,
    double* rhs,
    GridView grid,
    int keep_order,
    double gamma,
    double small_rho,
    double small_p,
    double inverse_dx,
    double inverse_dy,
    double inverse_dz,
    std::size_t physical_count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= physical_count) {
    return;
  }

  const int i = static_cast<int>(linear % grid.nx) + grid.nghost;
  const std::size_t jk = linear / grid.nx;
  const int j = static_cast<int>(jk % grid.ny) + grid.nghost;
  const int k = static_cast<int>(jk / grid.ny) + grid.nghost;
  const std::size_t center = cell_index(grid, i, j, k);

  constexpr double central6_coefficient[3] = {
      3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0};
  const int maximum_separation = keep_order == 2 ? 1 : 3;
  const double inverse_spacing[3] = {inverse_dx, inverse_dy, inverse_dz};
  double minus_flux[5];
  double plus_flux[5];
  double result[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

  for (int direction = 0; direction < 3; ++direction) {
    for (int separation = 1; separation <= maximum_separation; ++separation) {
      const int di = direction == 0 ? separation : 0;
      const int dj = direction == 1 ? separation : 0;
      const int dk = direction == 2 ? separation : 0;
      keep_flux(
          q, grid, cell_index(grid, i - di, j - dj, k - dk), center,
          direction, gamma, small_rho, small_p, minus_flux);
      keep_flux(
          q, grid, center, cell_index(grid, i + di, j + dj, k + dk),
          direction, gamma, small_rho, small_p, plus_flux);
      const double derivative_coefficient =
          keep_order == 2 ? 0.5
                          : central6_coefficient[separation - 1];
      const double weight =
          2.0 * derivative_coefficient * inverse_spacing[direction];
      for (int variable = 0; variable < 5; ++variable) {
        result[variable] -=
            weight * (plus_flux[variable] - minus_flux[variable]);
      }
    }
  }

  for (int variable = 0; variable < 5; ++variable) {
    rhs[center + grid.cell_count * static_cast<std::size_t>(variable)] =
        result[variable];
  }
}

#include "nse_cuda_weno5z_roe.cuh"
#include "nse_cuda_hybrid.cuh"

__global__ void viscous_central6_kernel(
    const double* primitive,
    double* rhs,
    GridView grid,
    double gamma,
    double inverse_reynolds,
    double inverse_prandtl,
    double inverse_dx,
    double inverse_dy,
    double inverse_dz,
    std::size_t physical_count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= physical_count) {
    return;
  }

  const int i = static_cast<int>(linear % grid.nx) + grid.nghost;
  const std::size_t jk = linear / grid.nx;
  const int j = static_cast<int>(jk % grid.ny) + grid.nghost;
  const int k = static_cast<int>(jk / grid.ny) + grid.nghost;
  const std::size_t center = cell_index(grid, i, j, k);
  const double* u_field = primitive;
  const double* v_field = primitive + grid.cell_count;
  const double* w_field = primitive + 2 * grid.cell_count;
  const double* temperature_field = primitive + 3 * grid.cell_count;
  const double inverse_dx2 = inverse_dx * inverse_dx;
  const double inverse_dy2 = inverse_dy * inverse_dy;
  const double inverse_dz2 = inverse_dz * inverse_dz;

  const double u = u_field[center];
  const double v = v_field[center];
  const double w = w_field[center];
  const double ux =
      first_derivative(u_field, grid, i, j, k, 0, inverse_dx);
  const double uy =
      first_derivative(u_field, grid, i, j, k, 1, inverse_dy);
  const double uz =
      first_derivative(u_field, grid, i, j, k, 2, inverse_dz);
  const double vx =
      first_derivative(v_field, grid, i, j, k, 0, inverse_dx);
  const double vy =
      first_derivative(v_field, grid, i, j, k, 1, inverse_dy);
  const double vz =
      first_derivative(v_field, grid, i, j, k, 2, inverse_dz);
  const double wx =
      first_derivative(w_field, grid, i, j, k, 0, inverse_dx);
  const double wy =
      first_derivative(w_field, grid, i, j, k, 1, inverse_dy);
  const double wz =
      first_derivative(w_field, grid, i, j, k, 2, inverse_dz);

  const double uxx =
      second_derivative(u_field, grid, i, j, k, 0, inverse_dx2);
  const double uyy =
      second_derivative(u_field, grid, i, j, k, 1, inverse_dy2);
  const double uzz =
      second_derivative(u_field, grid, i, j, k, 2, inverse_dz2);
  const double vxx =
      second_derivative(v_field, grid, i, j, k, 0, inverse_dx2);
  const double vyy =
      second_derivative(v_field, grid, i, j, k, 1, inverse_dy2);
  const double vzz =
      second_derivative(v_field, grid, i, j, k, 2, inverse_dz2);
  const double wxx =
      second_derivative(w_field, grid, i, j, k, 0, inverse_dx2);
  const double wyy =
      second_derivative(w_field, grid, i, j, k, 1, inverse_dy2);
  const double wzz =
      second_derivative(w_field, grid, i, j, k, 2, inverse_dz2);

  const double uxy = mixed_derivative(
      u_field, grid, i, j, k, 0, 1, inverse_dx * inverse_dy);
  const double uxz = mixed_derivative(
      u_field, grid, i, j, k, 0, 2, inverse_dx * inverse_dz);
  const double vxy = mixed_derivative(
      v_field, grid, i, j, k, 0, 1, inverse_dx * inverse_dy);
  const double vyz = mixed_derivative(
      v_field, grid, i, j, k, 1, 2, inverse_dy * inverse_dz);
  const double wxz = mixed_derivative(
      w_field, grid, i, j, k, 0, 2, inverse_dx * inverse_dz);
  const double wyz = mixed_derivative(
      w_field, grid, i, j, k, 1, 2, inverse_dy * inverse_dz);

  const double momentum_x =
      (4.0 / 3.0) * uxx + uyy + uzz + (vxy + wxz) / 3.0;
  const double momentum_y =
      vxx + (4.0 / 3.0) * vyy + vzz + (uxy + wyz) / 3.0;
  const double momentum_z =
      wxx + wyy + (4.0 / 3.0) * wzz + (uxz + vyz) / 3.0;

  const double div_velocity = ux + vy + wz;
  const double tau_xx = 2.0 * ux - (2.0 / 3.0) * div_velocity;
  const double tau_yy = 2.0 * vy - (2.0 / 3.0) * div_velocity;
  const double tau_zz = 2.0 * wz - (2.0 / 3.0) * div_velocity;
  const double tau_xy = uy + vx;
  const double tau_xz = uz + wx;
  const double tau_yz = vz + wy;
  const double dissipation =
      tau_xx * ux + tau_yy * vy + tau_zz * wz
      + tau_xy * (uy + vx) + tau_xz * (uz + wx)
      + tau_yz * (vz + wy);

  const double lap_temperature =
      second_derivative(
          temperature_field, grid, i, j, k, 0, inverse_dx2)
      + second_derivative(
          temperature_field, grid, i, j, k, 1, inverse_dy2)
      + second_derivative(
          temperature_field, grid, i, j, k, 2, inverse_dz2);
  const double heat_coefficient =
      gamma * inverse_prandtl / (gamma - 1.0);

  rhs[center + grid.cell_count] += inverse_reynolds * momentum_x;
  rhs[center + 2 * grid.cell_count] += inverse_reynolds * momentum_y;
  rhs[center + 3 * grid.cell_count] += inverse_reynolds * momentum_z;
  rhs[center + 4 * grid.cell_count] += inverse_reynolds * (
      u * momentum_x + v * momentum_y + w * momentum_z
      + dissipation + heat_coefficient * lap_temperature);
}

#if defined(NSE_FORCING_CUFFT) || defined(NSE_FORCING_CUFFTMP)
__device__ inline double complex_norm_squared(cufftDoubleComplex value) {
  return value.x * value.x + value.y * value.y;
}

__device__ inline double conservative_velocity(
    const double* q,
    const GridView& grid,
    int i,
    int j,
    int k,
    int variable,
    double small_rho) {
  const std::size_t cell = cell_index(grid, i, j, k);
  return q[cell + grid.cell_count * static_cast<std::size_t>(variable)]
      / fmax(q[cell], small_rho);
}

__device__ inline double velocity_derivative(
    const double* q,
    const GridView& grid,
    int i,
    int j,
    int k,
    int variable,
    int direction,
    double inverse_spacing,
    double small_rho) {
  const int di = direction == 0 ? 1 : 0;
  const int dj = direction == 1 ? 1 : 0;
  const int dk = direction == 2 ? 1 : 0;
  return inverse_spacing / 60.0 * (
      -conservative_velocity(
          q, grid, i - 3 * di, j - 3 * dj, k - 3 * dk,
          variable, small_rho)
      + 9.0 * conservative_velocity(
          q, grid, i - 2 * di, j - 2 * dj, k - 2 * dk,
          variable, small_rho)
      - 45.0 * conservative_velocity(
          q, grid, i - di, j - dj, k - dk, variable, small_rho)
      + 45.0 * conservative_velocity(
          q, grid, i + di, j + dj, k + dk, variable, small_rho)
      - 9.0 * conservative_velocity(
          q, grid, i + 2 * di, j + 2 * dj, k + 2 * dk,
          variable, small_rho)
      + conservative_velocity(
          q, grid, i + 3 * di, j + 3 * dj, k + 3 * dk,
          variable, small_rho));
}

__global__ void forcing_weighted_velocity_kernel(
    const double* q,
    cufftDoubleComplex* spectral,
    GridView grid,
    double small_rho,
    std::size_t physical_count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= physical_count) {
    return;
  }
  const int x = static_cast<int>(linear % grid.nx) + grid.nghost;
  const std::size_t yz = linear / grid.nx;
  const int y = static_cast<int>(yz % grid.ny) + grid.nghost;
  const int z = static_cast<int>(yz / grid.ny) + grid.nghost;
  const std::size_t cell = cell_index(grid, x, y, z);
  const double inverse_sqrt_rho = 1.0 / sqrt(fmax(q[cell], small_rho));
  for (int component = 0; component < 3; ++component) {
    cufftDoubleComplex value;
    value.x = q[cell + grid.cell_count * (component + 1)]
        * inverse_sqrt_rho;
    value.y = 0.0;
    spectral[linear + physical_count * component] = value;
  }
}

__global__ void forcing_helmholtz_kernel(
    cufftDoubleComplex* spectral,
    cufftDoubleComplex* phi,
    double* energy_s,
    double* energy_d,
    int nx,
    int ny,
    int nz,
    double lx,
    double ly,
    double lz,
    int spectrum,
    double k_cutoff,
    std::size_t count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= count) {
    return;
  }
  const int ix = static_cast<int>(linear % nx);
  const std::size_t yz = linear / nx;
  const int iy = static_cast<int>(yz % ny);
  const int iz = static_cast<int>(yz / ny);
  const int mx = ix <= nx / 2 ? ix : ix - nx;
  const int my = iy <= ny / 2 ? iy : iy - ny;
  const int mz = iz <= nz / 2 ? iz : iz - nz;
  constexpr double two_pi = 6.283185307179586476925286766559;
  const double kx = two_pi * static_cast<double>(mx) / lx;
  const double ky = two_pi * static_cast<double>(my) / ly;
  const double kz = two_pi * static_cast<double>(mz) / lz;
  const double k2 = kx * kx + ky * ky + kz * kz;
  const bool retained = k2 > 0.0
      && (spectrum == 1 || sqrt(k2) < k_cutoff);

  const cufftDoubleComplex wx = spectral[linear];
  const cufftDoubleComplex wy = spectral[linear + count];
  const cufftDoubleComplex wz = spectral[linear + 2 * count];
  cufftDoubleComplex projection{0.0, 0.0};
  cufftDoubleComplex sx{0.0, 0.0};
  cufftDoubleComplex sy{0.0, 0.0};
  cufftDoubleComplex sz{0.0, 0.0};
  if (retained) {
    projection.x = (kx * wx.x + ky * wy.x + kz * wz.x) / k2;
    projection.y = (kx * wx.y + ky * wy.y + kz * wz.y) / k2;
    sx = {wx.x - kx * projection.x, wx.y - kx * projection.y};
    sy = {wy.x - ky * projection.x, wy.y - ky * projection.y};
    sz = {wz.x - kz * projection.x, wz.y - kz * projection.y};
  }
  phi[linear] = projection;
  spectral[linear] = sx;
  spectral[linear + count] = sy;
  spectral[linear + 2 * count] = sz;
  energy_s[linear] = complex_norm_squared(sx)
      + complex_norm_squared(sy) + complex_norm_squared(sz);
  energy_d[linear] = retained ? k2 * complex_norm_squared(projection) : 0.0;
}

__global__ void forcing_pressure_dilatation_kernel(
    const double* q,
    double* values,
    GridView grid,
    double gamma,
    double small_rho,
    double small_p,
    double inverse_dx,
    double inverse_dy,
    double inverse_dz,
    std::size_t physical_count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= physical_count) {
    return;
  }
  const int i = static_cast<int>(linear % grid.nx) + grid.nghost;
  const std::size_t jk = linear / grid.nx;
  const int j = static_cast<int>(jk % grid.ny) + grid.nghost;
  const int k = static_cast<int>(jk / grid.ny) + grid.nghost;
  const std::size_t cell = cell_index(grid, i, j, k);
  double rho;
  double u;
  double v;
  double w;
  double pressure;
  primitive_state(
      q, grid, cell, gamma, small_rho, small_p,
      rho, u, v, w, pressure);
  const double divergence =
      velocity_derivative(
          q, grid, i, j, k, 1, 0, inverse_dx, small_rho)
      + velocity_derivative(
          q, grid, i, j, k, 2, 1, inverse_dy, small_rho)
      + velocity_derivative(
          q, grid, i, j, k, 3, 2, inverse_dz, small_rho);
  values[linear] = pressure * divergence;
}

__global__ void forcing_add_kernel(
    const cufftDoubleComplex* physical,
    const double* q,
    double* rhs,
    GridView grid,
    double coefficient,
    double normalization,
    double small_rho,
    std::size_t physical_count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= physical_count) {
    return;
  }
  const int i = static_cast<int>(linear % grid.nx) + grid.nghost;
  const std::size_t jk = linear / grid.nx;
  const int j = static_cast<int>(jk % grid.ny) + grid.nghost;
  const int k = static_cast<int>(jk / grid.ny) + grid.nghost;
  const std::size_t cell = cell_index(grid, i, j, k);
  const double scale = coefficient * normalization
      * sqrt(fmax(q[cell], small_rho));
  for (int component = 0; component < 3; ++component) {
    rhs[cell + grid.cell_count * (component + 1)] += scale
        * physical[linear + physical_count * component].x;
  }
}

__global__ void forcing_build_dilatational_kernel(
    cufftDoubleComplex* spectral,
    const cufftDoubleComplex* phi,
    int nx,
    int ny,
    int nz,
    double lx,
    double ly,
    double lz,
    int spectrum,
    double k_cutoff,
    std::size_t count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= count) {
    return;
  }
  const int ix = static_cast<int>(linear % nx);
  const std::size_t yz = linear / nx;
  const int iy = static_cast<int>(yz % ny);
  const int iz = static_cast<int>(yz / ny);
  const int mx = ix <= nx / 2 ? ix : ix - nx;
  const int my = iy <= ny / 2 ? iy : iy - ny;
  const int mz = iz <= nz / 2 ? iz : iz - nz;
  constexpr double two_pi = 6.283185307179586476925286766559;
  const double kx = two_pi * static_cast<double>(mx) / lx;
  const double ky = two_pi * static_cast<double>(my) / ly;
  const double kz = two_pi * static_cast<double>(mz) / lz;
  const double k2 = kx * kx + ky * ky + kz * kz;
  const bool retained = k2 > 0.0
      && (spectrum == 1 || sqrt(k2) < k_cutoff);
  const cufftDoubleComplex projection = retained
      ? phi[linear] : cufftDoubleComplex{0.0, 0.0};
  spectral[linear] = {kx * projection.x, kx * projection.y};
  spectral[linear + count] = {ky * projection.x, ky * projection.y};
  spectral[linear + 2 * count] = {kz * projection.x, kz * projection.y};
}

#if defined(NSE_FORCING_CUFFTMP)
__global__ void forcing_weighted_velocity_component_kernel(
    const double* q,
    cufftDoubleComplex* physical,
    GridView grid,
    int component,
    double small_rho,
    std::size_t count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= count) return;
  const int x = static_cast<int>(linear % grid.nx) + grid.nghost;
  const std::size_t yz = linear / grid.nx;
  const int y = static_cast<int>(yz % grid.ny) + grid.nghost;
  const int z = static_cast<int>(yz / grid.ny) + grid.nghost;
  const std::size_t cell = cell_index(grid, x, y, z);
  const double inverse_sqrt_rho = 1.0 / sqrt(fmax(q[cell], small_rho));
  physical[linear] = {
      q[cell + grid.cell_count * static_cast<std::size_t>(component + 1)]
          * inverse_sqrt_rho,
      0.0};
}

__global__ void forcing_helmholtz_distributed_kernel(
    cufftDoubleComplex* spectral_x,
    cufftDoubleComplex* spectral_y,
    cufftDoubleComplex* spectral_z,
    cufftDoubleComplex* phi,
    double* energy_s,
    double* energy_d,
    int nx,
    int local_ny,
    int global_ny,
    int global_nz,
    int y_start,
    int z_start,
    double lx,
    double ly,
    double lz,
    int spectrum,
    double k_cutoff,
    std::size_t count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= count) return;
  const int ix = static_cast<int>(linear % nx);
  const std::size_t yz = linear / nx;
  const int iy = static_cast<int>(yz % local_ny) + y_start;
  const int iz = static_cast<int>(yz / local_ny) + z_start;
  const int mx = ix <= nx / 2 ? ix : ix - nx;
  const int my = iy <= global_ny / 2 ? iy : iy - global_ny;
  const int mz = iz <= global_nz / 2 ? iz : iz - global_nz;
  constexpr double two_pi = 6.283185307179586476925286766559;
  const double kx = two_pi * static_cast<double>(mx) / lx;
  const double ky = two_pi * static_cast<double>(my) / ly;
  const double kz = two_pi * static_cast<double>(mz) / lz;
  const double k2 = kx * kx + ky * ky + kz * kz;
  const bool retained = k2 > 0.0
      && (spectrum == 1 || sqrt(k2) < k_cutoff);
  const cufftDoubleComplex wx = spectral_x[linear];
  const cufftDoubleComplex wy = spectral_y[linear];
  const cufftDoubleComplex wz = spectral_z[linear];
  cufftDoubleComplex projection{0.0, 0.0};
  cufftDoubleComplex sx{0.0, 0.0};
  cufftDoubleComplex sy{0.0, 0.0};
  cufftDoubleComplex sz{0.0, 0.0};
  if (retained) {
    projection.x = (kx * wx.x + ky * wy.x + kz * wz.x) / k2;
    projection.y = (kx * wx.y + ky * wy.y + kz * wz.y) / k2;
    sx = {wx.x - kx * projection.x, wx.y - kx * projection.y};
    sy = {wy.x - ky * projection.x, wy.y - ky * projection.y};
    sz = {wz.x - kz * projection.x, wz.y - kz * projection.y};
  }
  phi[linear] = projection;
  spectral_x[linear] = sx;
  spectral_y[linear] = sy;
  spectral_z[linear] = sz;
  energy_s[linear] = complex_norm_squared(sx)
      + complex_norm_squared(sy) + complex_norm_squared(sz);
  energy_d[linear] = retained ? k2 * complex_norm_squared(projection) : 0.0;
}

__global__ void forcing_build_dilatational_distributed_kernel(
    cufftDoubleComplex* spectral_x,
    cufftDoubleComplex* spectral_y,
    cufftDoubleComplex* spectral_z,
    const cufftDoubleComplex* phi,
    int nx,
    int local_ny,
    int global_ny,
    int global_nz,
    int y_start,
    int z_start,
    double lx,
    double ly,
    double lz,
    int spectrum,
    double k_cutoff,
    std::size_t count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= count) return;
  const int ix = static_cast<int>(linear % nx);
  const std::size_t yz = linear / nx;
  const int iy = static_cast<int>(yz % local_ny) + y_start;
  const int iz = static_cast<int>(yz / local_ny) + z_start;
  const int mx = ix <= nx / 2 ? ix : ix - nx;
  const int my = iy <= global_ny / 2 ? iy : iy - global_ny;
  const int mz = iz <= global_nz / 2 ? iz : iz - global_nz;
  constexpr double two_pi = 6.283185307179586476925286766559;
  const double kx = two_pi * static_cast<double>(mx) / lx;
  const double ky = two_pi * static_cast<double>(my) / ly;
  const double kz = two_pi * static_cast<double>(mz) / lz;
  const double k2 = kx * kx + ky * ky + kz * kz;
  const bool retained = k2 > 0.0
      && (spectrum == 1 || sqrt(k2) < k_cutoff);
  const cufftDoubleComplex projection = retained
      ? phi[linear] : cufftDoubleComplex{0.0, 0.0};
  spectral_x[linear] = {kx * projection.x, kx * projection.y};
  spectral_y[linear] = {ky * projection.x, ky * projection.y};
  spectral_z[linear] = {kz * projection.x, kz * projection.y};
}

__global__ void forcing_add_component_kernel(
    const cufftDoubleComplex* physical,
    const double* q,
    double* rhs,
    GridView grid,
    int component,
    double coefficient,
    double normalization,
    double small_rho,
    std::size_t count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= count) return;
  const int i = static_cast<int>(linear % grid.nx) + grid.nghost;
  const std::size_t yz = linear / grid.nx;
  const int j = static_cast<int>(yz % grid.ny) + grid.nghost;
  const int k = static_cast<int>(yz / grid.ny) + grid.nghost;
  const std::size_t cell = cell_index(grid, i, j, k);
  rhs[cell + grid.cell_count * static_cast<std::size_t>(component + 1)] +=
      coefficient * normalization * sqrt(fmax(q[cell], small_rho))
      * physical[linear].x;
}
#endif
#endif

__global__ void ssprk_stage_kernel(
    double* q,
    const double* q0,
    const double* rhs,
    GridView grid,
    double dt,
    int stage,
    std::size_t physical_count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= physical_count) {
    return;
  }

  const int i = static_cast<int>(linear % grid.nx) + grid.nghost;
  const std::size_t jk = linear / grid.nx;
  const int j = static_cast<int>(jk % grid.ny) + grid.nghost;
  const int k = static_cast<int>(jk / grid.ny) + grid.nghost;
  const std::size_t cell = cell_index(grid, i, j, k);
  for (int variable = 0; variable < 5; ++variable) {
    const std::size_t index =
        cell + grid.cell_count * static_cast<std::size_t>(variable);
    const double euler = q[index] + dt * rhs[index];
    if (stage == 1) {
      q[index] = q0[index] + dt * rhs[index];
    } else if (stage == 2) {
      q[index] = 0.75 * q0[index] + 0.25 * euler;
    } else {
      q[index] = (1.0 / 3.0) * q0[index] + (2.0 / 3.0) * euler;
    }
  }
}

dim3 block_count(std::size_t count) {
  constexpr unsigned int threads = 256;
  return dim3(static_cast<unsigned int>((count + threads - 1) / threads));
}

#if defined(NSE_FORCING_CUFFT)
bool reduce_forcing_sum(
    NseCudaContext* context,
    const double* values,
    double* result,
    const char* operation) {
  if (!check_cuda(
          cub::DeviceReduce::Sum(
              context->reduce_storage,
              context->reduce_storage_bytes,
              values,
              context->forcing_sum,
              static_cast<int>(context->physical_count)),
          operation)) {
    return false;
  }
  return check_cuda(
      cudaMemcpy(
          result,
          context->forcing_sum,
          sizeof(double),
          cudaMemcpyDeviceToHost),
      "copy forcing reduction to host");
}

bool launch_forcing(NseCudaContext* context) {
  if (!context->forcing_enabled) {
    return true;
  }
  const double point_count = static_cast<double>(context->physical_count);
  const double normalization = 1.0 / point_count;
  const double lx = context->grid.nx * context->dx;
  const double ly = context->grid.ny * context->dy;
  const double lz = context->grid.nz * context->dz;

  forcing_weighted_velocity_kernel<<<
      block_count(context->physical_count), 256>>>(
      context->q,
      context->forcing_spectral,
      context->grid,
      context->small_rho,
      context->physical_count);
  if (!check_cuda(cudaGetLastError(), "prepare weighted forcing velocity")
      || !check_cufft(
          cufftExecZ2Z(
              context->forcing_plan,
              context->forcing_spectral,
              context->forcing_spectral,
              CUFFT_FORWARD),
          "forward forcing FFT")) {
    return false;
  }

  forcing_helmholtz_kernel<<<
      block_count(context->physical_count), 256>>>(
      context->forcing_spectral,
      context->forcing_phi,
      context->speed,
      context->forcing_energy_d,
      context->grid.nx,
      context->grid.ny,
      context->grid.nz,
      lx,
      ly,
      lz,
      context->forcing_spectrum,
      context->forcing_k_cutoff,
      context->physical_count);
  if (!check_cuda(cudaGetLastError(), "forcing Helmholtz projection")) {
    return false;
  }

  double denominator_s = 0.0;
  double denominator_d = 0.0;
  if (!reduce_forcing_sum(
          context, context->speed, &denominator_s,
          "reduce solenoidal forcing energy")
      || !reduce_forcing_sum(
          context, context->forcing_energy_d, &denominator_d,
          "reduce dilatational forcing energy")) {
    return false;
  }
  denominator_s /= point_count * point_count;
  denominator_d /= point_count * point_count;

  forcing_pressure_dilatation_kernel<<<
      block_count(context->physical_count), 256>>>(
      context->q,
      context->speed,
      context->grid,
      context->gamma,
      context->small_rho,
      context->small_p,
      1.0 / context->dx,
      1.0 / context->dy,
      1.0 / context->dz,
      context->physical_count);
  if (!check_cuda(cudaGetLastError(), "forcing pressure dilatation")) {
    return false;
  }
  double pressure_dilatation = 0.0;
  if (!reduce_forcing_sum(
          context, context->speed, &pressure_dilatation,
          "reduce pressure dilatation")) {
    return false;
  }
  pressure_dilatation /= point_count;

  const double ratio = context->forcing_dilatational_ratio;
  const double target_s =
      context->forcing_target_dissipation / (1.0 + ratio);
  const double target_d = context->forcing_target_dissipation - target_s;
  if (denominator_s <= context->forcing_denominator_floor) {
    set_error("solenoidal forcing denominator is too small");
    return false;
  }
  double coefficient_s = target_s / denominator_s;
  const double numerator_d = target_d - pressure_dilatation;
  double coefficient_d = 0.0;
  if (denominator_d <= context->forcing_denominator_floor) {
    if (target_d > context->forcing_denominator_floor) {
      set_error("dilatational forcing denominator is too small");
      return false;
    }
  } else {
    coefficient_d = numerator_d / denominator_d;
  }
  if (context->forcing_max_coefficient > 0.0) {
    coefficient_s = std::clamp(
        coefficient_s,
        -context->forcing_max_coefficient,
        context->forcing_max_coefficient);
    coefficient_d = std::clamp(
        coefficient_d,
        -context->forcing_max_coefficient,
        context->forcing_max_coefficient);
  }

  if (!check_cufft(
          cufftExecZ2Z(
              context->forcing_plan,
              context->forcing_spectral,
              context->forcing_spectral,
              CUFFT_INVERSE),
          "inverse solenoidal forcing FFT")) {
    return false;
  }
  forcing_add_kernel<<<block_count(context->physical_count), 256>>>(
      context->forcing_spectral,
      context->q,
      context->rhs,
      context->grid,
      coefficient_s,
      normalization,
      context->small_rho,
      context->physical_count);
  if (!check_cuda(cudaGetLastError(), "add solenoidal forcing")) {
    return false;
  }

  forcing_build_dilatational_kernel<<<
      block_count(context->physical_count), 256>>>(
      context->forcing_spectral,
      context->forcing_phi,
      context->grid.nx,
      context->grid.ny,
      context->grid.nz,
      lx,
      ly,
      lz,
      context->forcing_spectrum,
      context->forcing_k_cutoff,
      context->physical_count);
  if (!check_cuda(cudaGetLastError(), "build dilatational forcing spectrum")
      || !check_cufft(
          cufftExecZ2Z(
              context->forcing_plan,
              context->forcing_spectral,
              context->forcing_spectral,
              CUFFT_INVERSE),
          "inverse dilatational forcing FFT")) {
    return false;
  }
  forcing_add_kernel<<<block_count(context->physical_count), 256>>>(
      context->forcing_spectral,
      context->q,
      context->rhs,
      context->grid,
      coefficient_d,
      normalization,
      context->small_rho,
      context->physical_count);
  if (!check_cuda(cudaGetLastError(), "add dilatational forcing")) {
    return false;
  }

  ++context->forcing_evaluations;
  if (context->forcing_report_interval > 0
      && (context->forcing_evaluations == 1
          || context->forcing_evaluations
                  % context->forcing_report_interval == 0)) {
    std::printf(
        "# forcing %d %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
        context->forcing_evaluations,
        coefficient_s,
        coefficient_d,
        denominator_s,
        denominator_d,
        pressure_dilatation,
        target_s,
        target_d);
  }
  return true;
}
#endif

#if defined(NSE_FORCING_CUFFTMP)
bool reduce_forcing_sum_distributed(
    NseCudaContext* context,
    const double* values,
    double* result,
    const char* operation) {
  double local = 0.0;
  if (!check_cuda(
          cub::DeviceReduce::Sum(
              context->reduce_storage,
              context->reduce_storage_bytes,
              values,
              context->forcing_sum,
              static_cast<int>(context->physical_count)),
          operation)
      || !check_cuda(
          cudaMemcpy(
              &local,
              context->forcing_sum,
              sizeof(double),
              cudaMemcpyDeviceToHost),
          "copy local distributed forcing reduction")) {
    return false;
  }
  if (MPI_Allreduce(
          &local, result, 1, MPI_DOUBLE, MPI_SUM,
          context->forcing_comm) != MPI_SUCCESS) {
    set_error(std::string("MPI_Allreduce failed during ") + operation);
    return false;
  }
  return true;
}

bool launch_forcing_distributed(NseCudaContext* context) {
  if (!context->forcing_enabled) return true;
  if (!context->forcing_distributed_ready) {
    set_error("cuFFTMp forcing was not configured before time advancement");
    return false;
  }
  const double point_count = static_cast<double>(context->grid.nx)
      * static_cast<double>(context->global_ny)
      * static_cast<double>(context->global_nz);
  const double normalization = 1.0 / point_count;
  const double lx = context->grid.nx * context->dx;
  const double ly = context->global_ny * context->dy;
  const double lz = context->global_nz * context->dz;
  cufftDoubleComplex* physical[3] = {};
  cufftDoubleComplex* spectral[3] = {};

  for (int component = 0; component < 3; ++component) {
    physical[component] = descriptor_data(
        context->forcing_physical[component]);
    spectral[component] = descriptor_data(
        context->forcing_spectral[component]);
    if (physical[component] == nullptr || spectral[component] == nullptr) {
      return false;
    }
    forcing_weighted_velocity_component_kernel<<<
        block_count(context->physical_count), 256>>>(
        context->q,
        physical[component],
        context->grid,
        component,
        context->small_rho,
        context->physical_count);
    if (!check_cuda(
            cudaGetLastError(),
            "prepare distributed weighted forcing velocity")
        || !check_cufft(
            cufftXtExecDescriptor(
                context->forcing_plan,
                context->forcing_physical[component],
                context->forcing_spectral[component],
                CUFFT_FORWARD),
            "forward distributed forcing FFT")) {
      return false;
    }
  }

  forcing_helmholtz_distributed_kernel<<<
      block_count(context->physical_count), 256>>>(
      spectral[0], spectral[1], spectral[2],
      context->forcing_phi,
      context->speed,
      context->forcing_energy_d,
      context->grid.nx,
      context->grid.ny,
      context->global_ny,
      context->global_nz,
      context->global_y_start,
      context->global_z_start,
      lx, ly, lz,
      context->forcing_spectrum,
      context->forcing_k_cutoff,
      context->physical_count);
  if (!check_cuda(
          cudaGetLastError(), "distributed forcing Helmholtz projection")) {
    return false;
  }

  double denominator_s = 0.0;
  double denominator_d = 0.0;
  if (!reduce_forcing_sum_distributed(
          context, context->speed, &denominator_s,
          "reduce distributed solenoidal forcing energy")
      || !reduce_forcing_sum_distributed(
          context, context->forcing_energy_d, &denominator_d,
          "reduce distributed dilatational forcing energy")) {
    return false;
  }
  denominator_s /= point_count * point_count;
  denominator_d /= point_count * point_count;

  forcing_pressure_dilatation_kernel<<<
      block_count(context->physical_count), 256>>>(
      context->q,
      context->speed,
      context->grid,
      context->gamma,
      context->small_rho,
      context->small_p,
      1.0 / context->dx,
      1.0 / context->dy,
      1.0 / context->dz,
      context->physical_count);
  if (!check_cuda(
          cudaGetLastError(), "distributed forcing pressure dilatation")) {
    return false;
  }
  double pressure_dilatation = 0.0;
  if (!reduce_forcing_sum_distributed(
          context, context->speed, &pressure_dilatation,
          "reduce distributed pressure dilatation")) {
    return false;
  }
  pressure_dilatation /= point_count;

  const double ratio = context->forcing_dilatational_ratio;
  const double target_s =
      context->forcing_target_dissipation / (1.0 + ratio);
  const double target_d = context->forcing_target_dissipation - target_s;
  if (denominator_s <= context->forcing_denominator_floor) {
    set_error("distributed solenoidal forcing denominator is too small");
    return false;
  }
  double coefficient_s = target_s / denominator_s;
  double coefficient_d = 0.0;
  if (denominator_d <= context->forcing_denominator_floor) {
    if (target_d > context->forcing_denominator_floor) {
      set_error("distributed dilatational forcing denominator is too small");
      return false;
    }
  } else {
    coefficient_d = (target_d - pressure_dilatation) / denominator_d;
  }
  if (context->forcing_max_coefficient > 0.0) {
    coefficient_s = std::clamp(
        coefficient_s,
        -context->forcing_max_coefficient,
        context->forcing_max_coefficient);
    coefficient_d = std::clamp(
        coefficient_d,
        -context->forcing_max_coefficient,
        context->forcing_max_coefficient);
  }

  for (int component = 0; component < 3; ++component) {
    if (!check_cufft(
            cufftXtExecDescriptor(
                context->forcing_plan,
                context->forcing_spectral[component],
                context->forcing_physical[component],
                CUFFT_INVERSE),
            "inverse distributed solenoidal forcing FFT")) {
      return false;
    }
    forcing_add_component_kernel<<<
        block_count(context->physical_count), 256>>>(
        physical[component], context->q, context->rhs, context->grid,
        component, coefficient_s, normalization, context->small_rho,
        context->physical_count);
    if (!check_cuda(
            cudaGetLastError(), "add distributed solenoidal forcing")) {
      return false;
    }
  }

  forcing_build_dilatational_distributed_kernel<<<
      block_count(context->physical_count), 256>>>(
      spectral[0], spectral[1], spectral[2], context->forcing_phi,
      context->grid.nx, context->grid.ny,
      context->global_ny, context->global_nz,
      context->global_y_start, context->global_z_start,
      lx, ly, lz,
      context->forcing_spectrum, context->forcing_k_cutoff,
      context->physical_count);
  if (!check_cuda(
          cudaGetLastError(),
          "build distributed dilatational forcing spectrum")) {
    return false;
  }
  for (int component = 0; component < 3; ++component) {
    if (!check_cufft(
            cufftXtExecDescriptor(
                context->forcing_plan,
                context->forcing_spectral[component],
                context->forcing_physical[component],
                CUFFT_INVERSE),
            "inverse distributed dilatational forcing FFT")) {
      return false;
    }
    forcing_add_component_kernel<<<
        block_count(context->physical_count), 256>>>(
        physical[component], context->q, context->rhs, context->grid,
        component, coefficient_d, normalization, context->small_rho,
        context->physical_count);
    if (!check_cuda(
            cudaGetLastError(), "add distributed dilatational forcing")) {
      return false;
    }
  }

  ++context->forcing_evaluations;
  if (context->forcing_rank == 0 && context->forcing_report_interval > 0
      && (context->forcing_evaluations == 1
          || context->forcing_evaluations
                  % context->forcing_report_interval == 0)) {
    std::printf(
        "# forcing %d %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
        context->forcing_evaluations,
        coefficient_s, coefficient_d,
        denominator_s, denominator_d,
        pressure_dilatation, target_s, target_d);
  }
  return true;
}
#endif

std::size_t halo_value_count(
    const NseCudaContext* context, int direction) {
  if (direction == halo_direction_y) {
    return static_cast<std::size_t>(context->grid.nx_total)
        * static_cast<std::size_t>(context->grid.nghost)
        * static_cast<std::size_t>(context->grid.nz_total)
        * static_cast<std::size_t>(context->nvar);
  }
  if (direction == halo_direction_z) {
    return static_cast<std::size_t>(context->grid.nx_total)
        * static_cast<std::size_t>(context->grid.ny_total)
        * static_cast<std::size_t>(context->grid.nghost)
        * static_cast<std::size_t>(context->nvar);
  }
  return 0;
}

bool valid_halo_request(
    const NseCudaContext* context, int direction, int side) {
  if (context == nullptr
      || (direction != halo_direction_y && direction != halo_direction_z)
      || (side != halo_side_low && side != halo_side_high)) {
    set_error("invalid CUDA halo-exchange request");
    return false;
  }
  if ((direction == halo_direction_y && !context->distributed_y)
      || (direction == halo_direction_z && !context->distributed_z)) {
    set_error("CUDA halo exchange requested for a local periodic direction");
    return false;
  }
  return true;
}

bool launch_periodic(NseCudaContext* context) {
  periodic_halo_kernel<<<block_count(context->grid.cell_count), 256>>>(
      context->q, context->grid, context->nvar,
      context->distributed_y ? 0 : 1,
      context->distributed_z ? 0 : 1);
  return check_cuda(cudaGetLastError(), "periodic halo kernel");
}

bool launch_rhs(NseCudaContext* context) {
  if (context->convective_scheme == convective_hybrid) {
    rhs_hybrid_kernel<<<block_count(context->physical_count), 256>>>(
        context->q,
        context->rhs,
        context->grid,
        context->hybrid_smooth_scheme,
        context->hybrid_shock_scheme,
        context->hybrid_sensor,
        context->hybrid_sensor_onset,
        context->hybrid_sensor_full,
        context->gamma,
        context->small_rho,
        context->small_p,
        1.0 / context->dx,
        1.0 / context->dy,
        1.0 / context->dz,
        context->physical_count);
  } else if (context->convective_scheme == convective_weno5z_roe) {
    rhs_weno5z_roe_kernel<<<block_count(context->physical_count), 256>>>(
        context->q,
        context->rhs,
        context->grid,
        context->gamma,
        context->small_rho,
        context->small_p,
        1.0 / context->dx,
        1.0 / context->dy,
        1.0 / context->dz,
        context->physical_count);
  } else {
    const int keep_order =
        context->convective_scheme == convective_keep2 ? 2 : 6;
    rhs_keep_kernel<<<block_count(context->physical_count), 256>>>(
        context->q,
        context->rhs,
        context->grid,
        keep_order,
        context->gamma,
        context->small_rho,
        context->small_p,
        1.0 / context->dx,
        1.0 / context->dy,
        1.0 / context->dz,
        context->physical_count);
  }
  if (!check_cuda(cudaGetLastError(), "convective right-hand-side kernel")) {
    return false;
  }
  if (context->viscous_enabled) {
    primitive_kernel<<<block_count(context->grid.cell_count), 256>>>(
        context->q,
        context->primitive,
        context->grid,
        context->gamma,
        context->small_rho,
        context->small_p);
    if (!check_cuda(cudaGetLastError(), "primitive-state kernel")) {
      return false;
    }
    viscous_central6_kernel<<<block_count(context->physical_count), 256>>>(
        context->primitive,
        context->rhs,
        context->grid,
        context->gamma,
        1.0 / context->reynolds,
        1.0 / context->prandtl,
        1.0 / context->dx,
        1.0 / context->dy,
        1.0 / context->dz,
        context->physical_count);
    if (!check_cuda(cudaGetLastError(), "central6 viscous kernel")) {
      return false;
    }
  }
#if defined(NSE_FORCING_CUFFT)
  return launch_forcing(context);
#elif defined(NSE_FORCING_CUFFTMP)
  return launch_forcing_distributed(context);
#else
  if (context->forcing_enabled) {
    set_error("forcing requested, but this CUDA build has no FFT backend");
    return false;
  }
  return true;
#endif
}

bool launch_stage(NseCudaContext* context, double dt, int stage) {
  ssprk_stage_kernel<<<block_count(context->physical_count), 256>>>(
      context->q,
      context->q0,
      context->rhs,
      context->grid,
      dt,
      stage,
      context->physical_count);
  return check_cuda(cudaGetLastError(), "SSPRK3 stage kernel");
}

}  // namespace

#if defined(NSE_INIT_CUFFT)
int transform_cuda_hit_complex_3d(
    int nx,
    int ny,
    int nz,
    int device,
    cufftDoubleComplex* host_field,
    int direction) {
  last_error.clear();
  if (nx <= 0 || ny <= 0 || nz <= 0 || host_field == nullptr
      || (direction != CUFFT_FORWARD && direction != CUFFT_INVERSE)) {
    set_error("invalid cuFFT HIT initialization argument");
    return 1;
  }

  const bool inverse = direction == CUFFT_INVERSE;
  const char* plan_error = inverse
      ? "create CUDA HIT inverse FFT plan"
      : "create CUDA HIT forward FFT plan";
  const char* execution_error = inverse
      ? "execute CUDA HIT inverse FFT"
      : "execute CUDA HIT forward FFT";

  const std::size_t nx_size = static_cast<std::size_t>(nx);
  const std::size_t ny_size = static_cast<std::size_t>(ny);
  const std::size_t nz_size = static_cast<std::size_t>(nz);
  if (nx_size > std::numeric_limits<std::size_t>::max() / ny_size
      || nx_size * ny_size
          > std::numeric_limits<std::size_t>::max() / nz_size) {
    set_error("cuFFT HIT grid size overflows size_t");
    return 1;
  }
  const std::size_t count = nx_size * ny_size * nz_size;
  if (count > std::numeric_limits<std::size_t>::max()
          / sizeof(cufftDoubleComplex)) {
    set_error("cuFFT HIT buffer size overflows size_t");
    return 1;
  }
  const std::size_t bytes = count * sizeof(cufftDoubleComplex);

  cufftDoubleComplex* device_field = nullptr;
  cufftHandle plan = 0;
  const auto cleanup = [&]() {
    if (plan != 0) {
      cufftDestroy(plan);
    }
    cudaFree(device_field);
  };

  if (!check_cuda(cudaSetDevice(device), "select CUDA HIT device")
      || !check_cuda(
          cudaMalloc(reinterpret_cast<void**>(&device_field), bytes),
          "allocate CUDA HIT Fourier field")
      || !check_cuda(
          cudaMemcpy(
              device_field, host_field, bytes, cudaMemcpyHostToDevice),
          "upload CUDA HIT Fourier field")
      || !check_cufft(
          cufftPlan3d(&plan, nz, ny, nx, CUFFT_Z2Z),
          plan_error)
      || !check_cufft(
          cufftExecZ2Z(
              plan, device_field, device_field, direction),
          execution_error)
      || !check_cuda(
          cudaMemcpy(
              host_field, device_field, bytes, cudaMemcpyDeviceToHost),
          "download CUDA HIT velocity field")) {
    cleanup();
    return 1;
  }

  cleanup();
  return 0;
}

NSE_CUDA_EXPORT int nse_cuda_inverse_complex_3d(
    int nx,
    int ny,
    int nz,
    int device,
    cufftDoubleComplex* host_field) {
  return transform_cuda_hit_complex_3d(
      nx, ny, nz, device, host_field, CUFFT_INVERSE);
}

NSE_CUDA_EXPORT int nse_cuda_forward_complex_3d(
    int nx,
    int ny,
    int nz,
    int device,
    cufftDoubleComplex* host_field) {
  return transform_cuda_hit_complex_3d(
      nx, ny, nz, device, host_field, CUFFT_FORWARD);
}
#endif

#if defined(NSE_FORCING_CUFFTMP) || defined(NSE_INIT_CUFFTMP)
bool make_cufftmp_pencil_plan(
    cufftHandle* plan,
    cudaLibXtDesc** physical,
    cudaLibXtDesc** spectral,
    int nx,
    int ny,
    int nz,
    int y_start,
    int y_count,
    int z_start,
    int z_count,
    MPI_Comm* comm) {
  int shape[3] = {nz, ny, nx};
  const long long lower[3] = {z_start, y_start, 0};
  const long long upper[3] = {
      z_start + z_count, y_start + y_count, nx};
  const long long strides[3] = {
      static_cast<long long>(y_count) * nx,
      static_cast<long long>(nx),
      1};
  std::size_t workspace = 0;
  return check_cufft(cufftCreate(plan), "create cuFFTMp plan")
      && check_cufft(
          cufftMpMakePlanDecomposition(
              *plan, 3, shape,
              lower, upper, strides,
              lower, upper, strides,
              CUFFT_Z2Z, comm, CUFFT_COMM_MPI, &workspace),
          "create cuFFTMp custom Y-Z pencil plan")
      && check_cufft(
          cufftXtMalloc(
              *plan, physical, CUFFT_XT_FORMAT_DISTRIBUTED_INPUT),
          "allocate cuFFTMp physical descriptor")
      && check_cufft(
          cufftXtMalloc(
              *plan, spectral, CUFFT_XT_FORMAT_DISTRIBUTED_OUTPUT),
          "allocate cuFFTMp spectral descriptor");
}
#endif

NSE_CUDA_EXPORT int nse_cuda_set_device(int device) {
  last_error.clear();
  if (device < 0) {
    set_error("CUDA device index must be non-negative");
    return 1;
  }
  return check_cuda(cudaSetDevice(device), "select CUDA device") ? 0 : 1;
}

#if defined(NSE_INIT_CUFFTMP)
NSE_CUDA_EXPORT void nse_cufftmp_fft_finalize() {
  if (hit_cufftmp.spectral != nullptr) {
    cufftXtFree(hit_cufftmp.spectral);
  }
  if (hit_cufftmp.physical != nullptr) {
    cufftXtFree(hit_cufftmp.physical);
  }
  if (hit_cufftmp.plan != 0) cufftDestroy(hit_cufftmp.plan);
  if (hit_cufftmp.comm != MPI_COMM_NULL) {
    MPI_Comm_free(&hit_cufftmp.comm);
  }
  hit_cufftmp = HitCufftMpContext{};
}

NSE_CUDA_EXPORT int nse_cufftmp_fft_initialize(
    int nx,
    int ny,
    int nz,
    int ylo,
    int yhi,
    int zlo,
    int zhi,
    int communicator) {
  last_error.clear();
  if (hit_cufftmp.plan != 0) {
    set_error("cuFFTMp HIT plan is already initialized");
    return 1;
  }
  if (nx <= 0 || ny <= 0 || nz <= 0
      || ylo < 1 || yhi < ylo || yhi > ny
      || zlo < 1 || zhi < zlo || zhi > nz) {
    set_error("invalid cuFFTMp HIT pencil bounds");
    return 1;
  }
  MPI_Comm source = MPI_Comm_f2c(static_cast<MPI_Fint>(communicator));
  if (source == MPI_COMM_NULL
      || MPI_Comm_dup(source, &hit_cufftmp.comm) != MPI_SUCCESS) {
    set_error("duplicate MPI communicator for cuFFTMp HIT initialization");
    nse_cufftmp_fft_finalize();
    return 1;
  }
  hit_cufftmp.local_count = static_cast<std::size_t>(nx)
      * static_cast<std::size_t>(yhi - ylo + 1)
      * static_cast<std::size_t>(zhi - zlo + 1);
  if (!make_cufftmp_pencil_plan(
          &hit_cufftmp.plan,
          &hit_cufftmp.physical,
          &hit_cufftmp.spectral,
          nx, ny, nz,
          ylo - 1, yhi - ylo + 1,
          zlo - 1, zhi - zlo + 1,
          &hit_cufftmp.comm)) {
    nse_cufftmp_fft_finalize();
    return 1;
  }
  return 0;
}

NSE_CUDA_EXPORT int nse_cufftmp_fft_execute(
    const cufftDoubleComplex* input,
    cufftDoubleComplex* output,
    int direction) {
  last_error.clear();
  if (hit_cufftmp.plan == 0 || input == nullptr || output == nullptr
      || (direction != CUFFT_FORWARD && direction != CUFFT_INVERSE)) {
    set_error("invalid cuFFTMp HIT transform request");
    return 1;
  }
  cudaLibXtDesc* input_descriptor = direction == CUFFT_FORWARD
      ? hit_cufftmp.physical : hit_cufftmp.spectral;
  cudaLibXtDesc* output_descriptor = direction == CUFFT_FORWARD
      ? hit_cufftmp.spectral : hit_cufftmp.physical;
  if (!check_cufft(
          cufftXtMemcpy(
              hit_cufftmp.plan,
              input_descriptor,
              const_cast<cufftDoubleComplex*>(input),
              CUFFT_COPY_HOST_TO_DEVICE),
          "upload local cuFFTMp HIT pencil")
      || !check_cufft(
          cufftXtExecDescriptor(
              hit_cufftmp.plan,
              input_descriptor,
              output_descriptor,
              direction),
          "execute cuFFTMp HIT transform")
      || !check_cufft(
          cufftXtMemcpy(
              hit_cufftmp.plan,
              output,
              output_descriptor,
              CUFFT_COPY_DEVICE_TO_HOST),
          "download local cuFFTMp HIT pencil")) {
    return 1;
  }
  return 0;
}
#endif

NSE_CUDA_EXPORT int nse_cuda_create(
    void** handle,
    int nx,
    int ny,
    int nz,
    int nghost,
    int nvar,
    int device,
    int distributed_y,
    int distributed_z,
    int convective_scheme,
    int hybrid_smooth_scheme,
    int hybrid_shock_scheme,
    int hybrid_sensor,
    int viscous_enabled,
    int forcing_enabled,
    int forcing_spectrum,
    int forcing_report_interval,
    double gamma,
    double cfl,
    double small_rho,
    double small_p,
    double reynolds,
    double prandtl,
    double dx,
    double dy,
    double dz,
    double forcing_k_cutoff,
    double forcing_target_dissipation,
    double forcing_dilatational_ratio,
    double forcing_denominator_floor,
    double forcing_max_coefficient,
    double hybrid_sensor_onset,
    double hybrid_sensor_full) {
  last_error.clear();
  if (handle == nullptr) {
    set_error("CUDA context output pointer is null");
    return 1;
  }
  *handle = nullptr;
  if (nx <= 0 || ny <= 0 || nz <= 0 || nghost < 3 || nvar != 5) {
    set_error("CUDA convective schemes require three ghosts and five variables");
    return 1;
  }
  if (convective_scheme != convective_keep2
      && convective_scheme != convective_keep6
      && convective_scheme != convective_weno5z_roe
      && convective_scheme != convective_hybrid) {
    set_error(
        "CUDA convective scheme must be KEEP2, KEEP6, WENO5Z_ROE, or HYBRID");
    return 1;
  }
  const bool valid_smooth_scheme =
      hybrid_smooth_scheme == convective_keep2
      || hybrid_smooth_scheme == convective_keep6
      || hybrid_smooth_scheme == convective_weno5z_roe;
  const bool valid_shock_scheme =
      hybrid_shock_scheme == convective_keep2
      || hybrid_shock_scheme == convective_keep6
      || hybrid_shock_scheme == convective_weno5z_roe;
  if (convective_scheme == convective_hybrid
      && (!valid_smooth_scheme || !valid_shock_scheme
          || hybrid_sensor != sensor_ducros_pressure
          || hybrid_sensor_onset < 0.0
          || hybrid_sensor_full <= hybrid_sensor_onset)) {
    set_error("invalid CUDA hybrid convective configuration");
    return 1;
  }
  if (gamma <= 1.0 || cfl <= 0.0 || dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
    set_error("invalid NSE physical or grid parameter");
    return 1;
  }
  if (viscous_enabled != 0
      && (nghost < 3 || reynolds <= 0.0 || prandtl <= 0.0)) {
    set_error("central6 viscosity requires three ghosts and positive Re/Pr");
    return 1;
  }
  if (forcing_enabled != 0) {
#if !defined(NSE_FORCING_CUFFTMP)
    if (distributed_y != 0 || distributed_z != 0) {
      set_error(
          "MPI+CUDA forcing requires a distributed FFT backend; "
          "the cuFFT backend is single-GPU only");
      return 1;
    }
#endif
#if !defined(NSE_FORCING_CUFFT) && !defined(NSE_FORCING_CUFFTMP)
    set_error("forcing requested, but this CUDA build has no FFT backend");
    return 1;
#else
    if ((forcing_spectrum != 1 && forcing_spectrum != 2)
        || forcing_target_dissipation <= 0.0
        || forcing_dilatational_ratio < 0.0
        || forcing_denominator_floor <= 0.0
        || forcing_max_coefficient < 0.0
        || forcing_report_interval < 0
        || (forcing_spectrum == 2 && forcing_k_cutoff <= 0.0)) {
      set_error("invalid Petersen-Livescu forcing parameter");
      return 1;
    }
#endif
  }

  NseCudaContext* context = new (std::nothrow) NseCudaContext();
  if (context == nullptr) {
    set_error("failed to allocate the CUDA context on the host");
    return 1;
  }
  context->device = device;
  context->distributed_y = distributed_y != 0;
  context->distributed_z = distributed_z != 0;
  context->global_ny = ny;
  context->global_nz = nz;
  context->global_y_start = 0;
  context->global_z_start = 0;
  context->grid.nx = nx;
  context->grid.ny = ny;
  context->grid.nz = nz;
  context->grid.nghost = nghost;
  context->grid.nx_total = nx + 2 * nghost;
  context->grid.ny_total = ny + 2 * nghost;
  context->grid.nz_total = nz + 2 * nghost;
  context->grid.cell_count =
      static_cast<std::size_t>(context->grid.nx_total)
      * static_cast<std::size_t>(context->grid.ny_total)
      * static_cast<std::size_t>(context->grid.nz_total);
  context->physical_count =
      static_cast<std::size_t>(nx)
      * static_cast<std::size_t>(ny)
      * static_cast<std::size_t>(nz);
  if (context->physical_count >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    set_error("CUDA reduction currently supports at most INT_MAX cells");
    release_context(context);
    return 1;
  }
  context->nvar = nvar;
  context->convective_scheme = convective_scheme;
  context->hybrid_smooth_scheme = hybrid_smooth_scheme;
  context->hybrid_shock_scheme = hybrid_shock_scheme;
  context->hybrid_sensor = hybrid_sensor;
  context->hybrid_sensor_onset = hybrid_sensor_onset;
  context->hybrid_sensor_full = hybrid_sensor_full;
  context->gamma = gamma;
  context->cfl = cfl;
  context->small_rho = small_rho;
  context->small_p = small_p;
  context->reynolds = reynolds;
  context->prandtl = prandtl;
  context->viscous_enabled = viscous_enabled != 0;
  context->forcing_enabled = forcing_enabled != 0;
  context->forcing_spectrum = forcing_spectrum;
  context->forcing_report_interval = forcing_report_interval;
  context->forcing_k_cutoff = forcing_k_cutoff;
  context->forcing_target_dissipation = forcing_target_dissipation;
  context->forcing_dilatational_ratio = forcing_dilatational_ratio;
  context->forcing_denominator_floor = forcing_denominator_floor;
  context->forcing_max_coefficient = forcing_max_coefficient;
  context->dx = dx;
  context->dy = dy;
  context->dz = dz;
  context->state_bytes =
      context->grid.cell_count * static_cast<std::size_t>(nvar)
      * sizeof(double);
  const std::size_t y_halo_count = halo_value_count(
      context, halo_direction_y);
  const std::size_t z_halo_count = halo_value_count(
      context, halo_direction_z);
  context->halo_buffer_count = std::max(
      context->distributed_y ? y_halo_count : std::size_t{0},
      context->distributed_z ? z_halo_count : std::size_t{0});

  if (!check_cuda(cudaSetDevice(device), "select CUDA device")
      || !check_cuda(cudaFree(nullptr), "initialize CUDA runtime")
      || !check_cuda(
          cudaMalloc(reinterpret_cast<void**>(&context->q),
                     context->state_bytes),
          "allocate Q")
      || !check_cuda(
          cudaMalloc(reinterpret_cast<void**>(&context->q0),
                     context->state_bytes),
          "allocate Q0")
      || !check_cuda(
          cudaMalloc(reinterpret_cast<void**>(&context->rhs),
                     context->state_bytes),
          "allocate RHS")
      || (context->viscous_enabled
          && !check_cuda(
              cudaMalloc(reinterpret_cast<void**>(&context->primitive),
                         4 * context->grid.cell_count * sizeof(double)),
              "allocate primitive workspace"))
      || !check_cuda(
          cudaMalloc(reinterpret_cast<void**>(&context->speed),
                     context->physical_count * sizeof(double)),
          "allocate CFL speed")
      || !check_cuda(
          cudaMalloc(reinterpret_cast<void**>(&context->max_speed),
                     sizeof(double)),
          "allocate maximum CFL speed")
      || (context->halo_buffer_count > 0
          && !check_cuda(
              cudaMalloc(reinterpret_cast<void**>(&context->halo_buffer),
                         context->halo_buffer_count * sizeof(double)),
              "allocate MPI halo staging buffer"))) {
    release_context(context);
    return 1;
  }

#if defined(NSE_FORCING_CUFFT)
  if (context->forcing_enabled) {
    const std::size_t spectral_bytes = 3 * context->physical_count
        * sizeof(cufftDoubleComplex);
    const std::size_t scalar_spectral_bytes = context->physical_count
        * sizeof(cufftDoubleComplex);
    if (!check_cuda(
            cudaMalloc(
                reinterpret_cast<void**>(&context->forcing_spectral),
                spectral_bytes),
            "allocate forcing spectra")
        || !check_cuda(
            cudaMalloc(
                reinterpret_cast<void**>(&context->forcing_phi),
                scalar_spectral_bytes),
            "allocate forcing Helmholtz potential")
        || !check_cuda(
            cudaMalloc(
                reinterpret_cast<void**>(&context->forcing_energy_d),
                context->physical_count * sizeof(double)),
            "allocate dilatational forcing energy")
        || !check_cuda(
            cudaMalloc(
                reinterpret_cast<void**>(&context->forcing_sum),
                sizeof(double)),
            "allocate forcing reduction result")) {
      release_context(context);
      return 1;
    }
    int dimensions[3] = {nz, ny, nx};
    const int distance = static_cast<int>(context->physical_count);
    if (!check_cufft(
            cufftPlanMany(
                &context->forcing_plan,
                3,
                dimensions,
                nullptr,
                1,
                distance,
                nullptr,
                1,
                distance,
                CUFFT_Z2Z,
                3),
            "create batched forcing FFT plan")) {
      release_context(context);
      return 1;
    }
  }
#endif

  std::size_t max_reduce_bytes = 0;
  cudaError_t status = cub::DeviceReduce::Max(
      nullptr,
      max_reduce_bytes,
      context->speed,
      context->max_speed,
      static_cast<int>(context->physical_count));
  if (!check_cuda(status, "query CUB maximum reduction workspace")) {
    release_context(context);
    return 1;
  }
  context->reduce_storage_bytes = max_reduce_bytes;
#if defined(NSE_FORCING_CUFFT)
  if (context->forcing_enabled) {
    std::size_t sum_reduce_bytes = 0;
    status = cub::DeviceReduce::Sum(
        nullptr,
        sum_reduce_bytes,
        context->speed,
        context->forcing_sum,
        static_cast<int>(context->physical_count));
    if (!check_cuda(status, "query CUB sum reduction workspace")) {
      release_context(context);
      return 1;
    }
    context->reduce_storage_bytes =
        std::max(context->reduce_storage_bytes, sum_reduce_bytes);
  }
#endif
  if (!check_cuda(
          cudaMalloc(&context->reduce_storage, context->reduce_storage_bytes),
          "allocate CUB reduction workspace")) {
    release_context(context);
    return 1;
  }

  *handle = context;
  return 0;
}

NSE_CUDA_EXPORT int nse_cuda_configure_cufftmp(
    void* handle,
    int global_ny,
    int global_nz,
    int global_y_start,
    int global_z_start,
    int communicator) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (context == nullptr) {
    set_error("invalid CUDA context during cuFFTMp configuration");
    return 1;
  }
#if !defined(NSE_FORCING_CUFFTMP)
  static_cast<void>(global_ny);
  static_cast<void>(global_nz);
  static_cast<void>(global_y_start);
  static_cast<void>(global_z_start);
  static_cast<void>(communicator);
  set_error("this CUDA build does not include the cuFFTMp forcing backend");
  return 1;
#else
  if (!context->forcing_enabled) {
    set_error("cuFFTMp configuration requested while forcing is disabled");
    return 1;
  }
  if (context->forcing_distributed_ready) {
    set_error("cuFFTMp forcing is already configured");
    return 1;
  }
  if (global_ny <= 0 || global_nz <= 0
      || global_y_start < 0 || global_z_start < 0
      || global_y_start + context->grid.ny > global_ny
      || global_z_start + context->grid.nz > global_nz) {
    set_error("invalid global Y-Z pencil bounds for cuFFTMp forcing");
    return 1;
  }
  MPI_Comm source = MPI_Comm_f2c(static_cast<MPI_Fint>(communicator));
  if (source == MPI_COMM_NULL
      || MPI_Comm_dup(source, &context->forcing_comm) != MPI_SUCCESS
      || MPI_Comm_rank(
          context->forcing_comm, &context->forcing_rank) != MPI_SUCCESS) {
    set_error("initialize MPI communicator for cuFFTMp forcing");
    return 1;
  }
  context->global_ny = global_ny;
  context->global_nz = global_nz;
  context->global_y_start = global_y_start;
  context->global_z_start = global_z_start;

  if (!make_cufftmp_pencil_plan(
          &context->forcing_plan,
          &context->forcing_physical[0],
          &context->forcing_spectral[0],
          context->grid.nx, global_ny, global_nz,
          global_y_start, context->grid.ny,
          global_z_start, context->grid.nz,
          &context->forcing_comm)) {
    return 1;
  }
  for (int component = 1; component < 3; ++component) {
    if (!check_cufft(
            cufftXtMalloc(
                context->forcing_plan,
                &context->forcing_physical[component],
                CUFFT_XT_FORMAT_DISTRIBUTED_INPUT),
            "allocate cuFFTMp forcing physical descriptor")
        || !check_cufft(
            cufftXtMalloc(
                context->forcing_plan,
                &context->forcing_spectral[component],
                CUFFT_XT_FORMAT_DISTRIBUTED_OUTPUT),
            "allocate cuFFTMp forcing spectral descriptor")) {
      return 1;
    }
  }
  if (!check_cuda(
          cudaMalloc(
              reinterpret_cast<void**>(&context->forcing_phi),
              context->physical_count * sizeof(cufftDoubleComplex)),
          "allocate distributed forcing Helmholtz potential")
      || !check_cuda(
          cudaMalloc(
              reinterpret_cast<void**>(&context->forcing_energy_d),
              context->physical_count * sizeof(double)),
          "allocate distributed forcing dilatational energy")
      || !check_cuda(
          cudaMalloc(
              reinterpret_cast<void**>(&context->forcing_sum),
              sizeof(double)),
          "allocate distributed forcing reduction result")) {
    return 1;
  }

  std::size_t sum_reduce_bytes = 0;
  if (!check_cuda(
          cub::DeviceReduce::Sum(
              nullptr,
              sum_reduce_bytes,
              context->speed,
              context->forcing_sum,
              static_cast<int>(context->physical_count)),
          "query distributed forcing reduction workspace")) {
    return 1;
  }
  if (sum_reduce_bytes > context->reduce_storage_bytes) {
    cudaFree(context->reduce_storage);
    context->reduce_storage = nullptr;
    context->reduce_storage_bytes = sum_reduce_bytes;
    if (!check_cuda(
            cudaMalloc(
                &context->reduce_storage,
                context->reduce_storage_bytes),
            "resize distributed forcing reduction workspace")) {
      return 1;
    }
  }
  context->forcing_distributed_ready = true;
  if (context->forcing_rank == 0) {
    std::printf("# Petersen-Livescu forcing initialized\n");
    std::printf("# forcing FFT backend: cuFFTMp custom Y-Z pencils\n");
  }
  return 0;
#endif
}

NSE_CUDA_EXPORT int nse_cuda_upload(void* handle, const double* q) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (context == nullptr || q == nullptr) {
    set_error("invalid CUDA context or upload pointer");
    return 1;
  }
  if (!check_cuda(
          cudaMemcpy(
              context->q, q, context->state_bytes, cudaMemcpyHostToDevice),
          "copy Q to CUDA device")
      || !launch_periodic(context)
      || !check_cuda(cudaDeviceSynchronize(), "finish NSE upload")) {
    return 1;
  }
  return 0;
}

NSE_CUDA_EXPORT int nse_cuda_download(void* handle, double* q) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (context == nullptr || q == nullptr) {
    set_error("invalid CUDA context or download pointer");
    return 1;
  }
  if (!launch_periodic(context)
      || !check_cuda(
          cudaMemcpy(
              q, context->q, context->state_bytes, cudaMemcpyDeviceToHost),
          "copy Q from CUDA device")) {
    return 1;
  }
  return 0;
}

NSE_CUDA_EXPORT int nse_cuda_apply_local_periodic(void* handle) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (context == nullptr) {
    set_error("invalid CUDA context during local periodic update");
    return 1;
  }
  return launch_periodic(context) ? 0 : 1;
}

NSE_CUDA_EXPORT int nse_cuda_halo_count(
    void* handle, int direction, std::size_t* count) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (context == nullptr || count == nullptr
      || (direction != halo_direction_y && direction != halo_direction_z)) {
    set_error("invalid CUDA halo-count request");
    return 1;
  }
  *count = halo_value_count(context, direction);
  return 0;
}

NSE_CUDA_EXPORT int nse_cuda_pack_halo(
    void* handle, int direction, int side, double* host_buffer) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (host_buffer == nullptr
      || !valid_halo_request(context, direction, side)) {
    if (host_buffer == nullptr) set_error("CUDA halo pack buffer is null");
    return 1;
  }
  const std::size_t count = halo_value_count(context, direction);
  if (count > context->halo_buffer_count || context->halo_buffer == nullptr) {
    set_error("CUDA halo staging buffer is too small");
    return 1;
  }
  if (direction == halo_direction_y) {
    pack_y_halo_kernel<<<block_count(count), 256>>>(
        context->q, context->halo_buffer, context->grid, context->nvar,
        side, count);
  } else {
    pack_z_halo_kernel<<<block_count(count), 256>>>(
        context->q, context->halo_buffer, context->grid, context->nvar,
        side, count);
  }
  if (!check_cuda(cudaGetLastError(), "pack CUDA MPI halo")
      || !check_cuda(
          cudaMemcpy(host_buffer, context->halo_buffer,
                     count * sizeof(double), cudaMemcpyDeviceToHost),
          "copy packed CUDA MPI halo to host")) {
    return 1;
  }
  return 0;
}

NSE_CUDA_EXPORT int nse_cuda_unpack_halo(
    void* handle, int direction, int side, const double* host_buffer) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (host_buffer == nullptr
      || !valid_halo_request(context, direction, side)) {
    if (host_buffer == nullptr) set_error("CUDA halo unpack buffer is null");
    return 1;
  }
  const std::size_t count = halo_value_count(context, direction);
  if (count > context->halo_buffer_count || context->halo_buffer == nullptr) {
    set_error("CUDA halo staging buffer is too small");
    return 1;
  }
  if (!check_cuda(
          cudaMemcpy(context->halo_buffer, host_buffer,
                     count * sizeof(double), cudaMemcpyHostToDevice),
          "copy received MPI halo to CUDA device")) {
    return 1;
  }
  if (direction == halo_direction_y) {
    unpack_y_halo_kernel<<<block_count(count), 256>>>(
        context->q, context->halo_buffer, context->grid, context->nvar,
        side, count);
  } else {
    unpack_z_halo_kernel<<<block_count(count), 256>>>(
        context->q, context->halo_buffer, context->grid, context->nvar,
        side, count);
  }
  return check_cuda(cudaGetLastError(), "unpack CUDA MPI halo") ? 0 : 1;
}

NSE_CUDA_EXPORT int nse_cuda_compute_dt(void* handle, double* dt) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (context == nullptr || dt == nullptr) {
    set_error("invalid CUDA context or time-step pointer");
    return 1;
  }

  wave_speed_kernel<<<block_count(context->physical_count), 256>>>(
      context->q,
      context->speed,
      context->grid,
      context->gamma,
      context->cfl,
      context->small_rho,
      context->small_p,
      context->reynolds,
      context->prandtl,
      1.0 / context->dx,
      1.0 / context->dy,
      1.0 / context->dz,
      context->viscous_enabled,
      context->physical_count);
  if (!check_cuda(cudaGetLastError(), "CFL wave-speed kernel")) {
    return 1;
  }
  if (!check_cuda(
          cub::DeviceReduce::Max(
              context->reduce_storage,
              context->reduce_storage_bytes,
              context->speed,
              context->max_speed,
              static_cast<int>(context->physical_count)),
          "reduce maximum CFL speed")) {
    return 1;
  }

  double max_speed = 0.0;
  if (!check_cuda(
          cudaMemcpy(
              &max_speed,
              context->max_speed,
              sizeof(double),
              cudaMemcpyDeviceToHost),
          "copy maximum CFL speed")) {
    return 1;
  }
  if (!std::isfinite(max_speed) || max_speed <= 0.0) {
    set_error("computed CFL wave speed is not positive and finite");
    return 1;
  }
  *dt = context->cfl / max_speed;
  return 0;
}

NSE_CUDA_EXPORT int nse_cuda_begin_ssprk3(void* handle, double dt) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (context == nullptr) {
    set_error("invalid CUDA context");
    return 1;
  }
  if (!std::isfinite(dt) || dt <= 0.0) {
    set_error("SSPRK3 time step must be positive and finite");
    return 1;
  }
  return check_cuda(
      cudaMemcpy(context->q0, context->q, context->state_bytes,
                 cudaMemcpyDeviceToDevice),
      "copy Q to Q0") ? 0 : 1;
}

NSE_CUDA_EXPORT int nse_cuda_advance_ssprk3_stage(
    void* handle, double dt, int stage) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (context == nullptr || !std::isfinite(dt) || dt <= 0.0
      || stage < 1 || stage > 3) {
    set_error("invalid CUDA SSPRK3 stage request");
    return 1;
  }
  if (!launch_rhs(context) || !launch_stage(context, dt, stage)) {
    return 1;
  }
  return 0;
}

NSE_CUDA_EXPORT int nse_cuda_advance_ssprk3(void* handle, double dt) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (context == nullptr) {
    set_error("invalid CUDA context");
    return 1;
  }
  if (nse_cuda_begin_ssprk3(handle, dt) != 0) {
    return 1;
  }

  for (int stage = 1; stage <= 3; ++stage) {
    if (!launch_periodic(context)
        || nse_cuda_advance_ssprk3_stage(handle, dt, stage) != 0) {
      return 1;
    }
  }
  return 0;
}

NSE_CUDA_EXPORT int nse_cuda_synchronize(void* handle) {
  last_error.clear();
  auto* context = static_cast<NseCudaContext*>(handle);
  if (context == nullptr) {
    set_error("invalid CUDA context");
    return 1;
  }
  return check_cuda(cudaDeviceSynchronize(), "synchronize CUDA device") ? 0 : 1;
}

NSE_CUDA_EXPORT void nse_cuda_destroy(void* handle) {
  release_context(static_cast<NseCudaContext*>(handle));
}

NSE_CUDA_EXPORT void nse_cuda_get_last_error(
    char* buffer, int buffer_size) {
  if (buffer == nullptr || buffer_size <= 0) {
    return;
  }
  const std::size_t count = std::min(
      last_error.size(), static_cast<std::size_t>(buffer_size - 1));
  std::memcpy(buffer, last_error.data(), count);
  buffer[count] = '\0';
}
