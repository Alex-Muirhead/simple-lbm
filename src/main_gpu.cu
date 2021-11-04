#include <array>
#include <iostream>
#include <sstream>

#include "hdfio.h"
#include "index.h"

using namespace std;

void checkError(cudaError_t e)
{
   if (e != cudaSuccess)
   {
      std::cerr << "CUDA error: " << int(e) << " : " << cudaGetErrorString(e) << '\n';
      abort();
   }
}

static const int D = 2; // 2 dimensions to the simulation
static const int Q = 9; // 9 different velocity components in the distribution

__device__ static const float c[Q][D] = {
    // Stationary (v* = 0)
    {  0,  0 },
    // Cardinal directions (v* = 1)
    {  0, -1 },
    {  1,  0 },
    {  0,  1 },
    { -1,  0 },
    // Diagonal directions (v* = sqrt{2})
    {  1, -1 },
    {  1,  1 },
    { -1,  1 },
    { -1, -1 }
};

__device__ static const float weights[Q] = {
    // Stationary
    4./9,
    // Cardinal directions
    1./9,
    1./9,
    1./9,
    1./9,
    // Diagonal directions
    1./36,
    1./36,
    1./36,
    1./36
};


enum PointType : int {
    Fluid = 1 << 0,   // Standard streaming
    NoSlip = 1 << 1,  // Bounce-back boundary
    Slip_V = 1 << 2,  // Vertical bounce-forward boundary
    Slip_H = 1 << 3,  // Horiztonal bounce-forward boundary
};

__device__ static unsigned int gpu_size_x;
__device__ static unsigned int gpu_size_y;

static unsigned int size_x;
static unsigned int size_y;
static unsigned int nThreads = 16;

__device__ __forceinline__
size_t scalar_index(int x, int y) {
    return gpu_size_x * y + x;
}

__device__ __forceinline__
size_t field_index(int x, int y, int q) {
    return Q * (gpu_size_x * y + x) + q;
}

__global__
void set_gpu_size(int size_x, int size_y) {
    gpu_size_x = size_x;
    gpu_size_y = size_y;
}

__global__
void init_kernel(double* f, double* rho, double* vel_x, double* vel_y) {
    unsigned int y = blockIdx.y;
    unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;

    if ((x >= gpu_size_x) || (y >= gpu_size_y)) {
        // Exit kernel if out of bounds
        return;
    }

    // Fetch macroscopic properties
    double u = vel_x[scalar_index(x, y)];
    double v = vel_y[scalar_index(x, y)];
    double r =   rho[scalar_index(x, y)];

    // Initialise equilibrium
    double vel_sqr = u * u + v * v;

    for (int q = 0; q < Q; q++) {
        double vel_dot = u * c[q][0] + v * c[q][1];
        f[field_index(x, y, q)] = weights[q] * r * (1. + 3.0 * vel_dot - 1.5 * vel_sqr + 4.5 * vel_dot * vel_dot);
    }
}

__global__
void collide_kernel(double* f, double omega) {
    unsigned int y = blockIdx.y;
    unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;

    if ((x >= gpu_size_x) || (y >= gpu_size_y)) {
        // Exit kernel if out of bounds
        return;
    }

    // -- Calculate macroscopic properties --

    double r = 0.0;  // Density
    double u = 0.0;  // Velocity x
    double v = 0.0;  // Velocity y

    for (int q = 0; q < Q; q++) {
        double dist = f[field_index(x, y, q)];

        r += dist;            // 0th order moment
        u += c[q][0] * dist;  // 1st order moment
        v += c[q][1] * dist;  // ^^^
    }
    // Velocity = Momentum / Density
    u /= r;
    v /= r;

    // Precompute velocity norm
    double vel_sqr = u * u + v * v;

    // Can unroll loop if necessary
    for (int q = 0; q < Q; q++) {
        double vel_dot = u * c[q][0] + v * c[q][1];
        double feq = weights[q] * r * (1. + 3.0 * vel_dot - 1.5 * vel_sqr + 4.5 * vel_dot * vel_dot);
        f[field_index(x, y, q)] = (1 - omega) * f[field_index(x, y, q)] + omega * feq;
    }
}

__global__
void stream_kernel(double* f_src, double* f_dst) {
    unsigned int y = blockIdx.y;
    unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;

    if ((x >= gpu_size_x) || (y >= gpu_size_y)) {
        // Exit kernel if out of bounds
        return;
    }

    // Standard streaming w/ periodic boundary
    for (unsigned int q = 0; q < Q; q++) {
        // Wrap coordinates over domain
        unsigned int xn = int(gpu_size_x + x + c[q][0]) % gpu_size_x;
        unsigned int yn = int(gpu_size_y + y + c[q][1]) % gpu_size_y;

        f_dst[field_index(xn, yn, q)] = f_src[field_index(x, y, q)];
    }
}

__global__
void properties_kernel(double* f, double* rho, double* vel_x, double* vel_y) {
    unsigned int y = blockIdx.y;
    unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;

    if ((x >= gpu_size_x) || (y >= gpu_size_y)) {
        // Exit kernel if out of bounds
        return;
    }

    // -- Calculate macroscopic properties --

    double r = 0.0;  // Density
    double u = 0.0;  // Velocity x
    double v = 0.0;  // Velocity y

    for (int q = 0; q < Q; q++) {
        double dist = f[field_index(x, y, q)];

        r += dist;            // 0th order moment
        u += c[q][0] * dist;  // 1st order moment
        v += c[q][1] * dist;  // ^^^
    }
    // Velocity = Momentum / Density
    u /= r;
    v /= r;

    rho  [scalar_index(x, y)] = r;
    vel_x[scalar_index(x, y)] = u;
    vel_y[scalar_index(x, y)] = v;
}

void init_eq(double* f, double* rho, double* vel_x, double* vel_y) {
    dim3 grid(size_x/nThreads, size_y, 1);  // blocks in grid
    dim3 threads(nThreads, 1, 1);           // threads in block

    init_kernel<<<grid, threads>>>(f, rho, vel_x, vel_y);
}

void collide(double* f, double omega) {
    dim3 grid(size_x/nThreads, size_y, 1);  // blocks in grid
    dim3 threads(nThreads, 1, 1);           // threads in block

    collide_kernel<<<grid, threads>>>(f, omega);
}

void stream(double* f_src, double* f_dst) {
    dim3 grid(size_x/nThreads, size_y, 1);  // blocks in grid
    dim3 threads(nThreads, 1, 1);           // threads in block

    stream_kernel<<<grid, threads>>>(f_src, f_dst);
}

void calculate_flow_properties(double* f, double* rho, double* vel_x, double* vel_y) {
    dim3 grid(size_x/nThreads, size_y, 1);  // blocks in grid
    dim3 threads(nThreads, 1, 1);           // threads in block

    properties_kernel<<<grid, threads>>>(f, rho, vel_x, vel_y);
}

int main(int argc, char* argv[]) {
    int val;
    if (argc >= 2) {
        std::istringstream iss(argv[1]);

        if (iss >> val) {
            // Conversion successful
        }
    }

    double nu_lb = 0.092;                     // Lattice dynamic viscosity
    double omega = 1.0 / (3. * nu_lb + 0.5);  // Relaxation parameter

    InputData input = InputData::open("taylor_green.h5");

    Field::set_field_shape(128, 128, Q);
    Scalar::set_scalar_shape(128, 128);

    size_x = 128;
    size_y = 128;

    set_gpu_size<<<1, 1>>>(size_x, size_y);

    int* type_lattice = new int[Scalar::size];

    input.configure_domain(type_lattice);

    double* rho   = new double[Scalar::size];
    double* vel_x = new double[Scalar::size];
    double* vel_y = new double[Scalar::size];

    double* rho_gpu;
    double* vel_x_gpu;
    double* vel_y_gpu;

    checkError(cudaMalloc((void**)&rho_gpu,   Scalar::size*sizeof(double)));
    checkError(cudaMalloc((void**)&vel_x_gpu, Scalar::size*sizeof(double)));
    checkError(cudaMalloc((void**)&vel_y_gpu, Scalar::size*sizeof(double)));

    input.read_scalar("initial/density", rho);
    input.read_scalar("initial/vel_x", vel_x);
    input.read_scalar("initial/vel_y", vel_y);

    checkError(cudaMemcpy(  rho_gpu,   rho, Scalar::size*sizeof(double), cudaMemcpyHostToDevice));
    checkError(cudaMemcpy(vel_x_gpu, vel_x, Scalar::size*sizeof(double), cudaMemcpyHostToDevice));
    checkError(cudaMemcpy(vel_y_gpu, vel_y, Scalar::size*sizeof(double), cudaMemcpyHostToDevice));

    input.close();

    double* f1_gpu;
    double* f2_gpu;

    checkError(cudaMalloc((void**)&f1_gpu, Field::size*sizeof(double)));
    checkError(cudaMalloc((void**)&f2_gpu, Field::size*sizeof(double)));

    // Set up equilibrium values for initial flow
    init_eq(f1_gpu, rho_gpu, vel_x_gpu, vel_y_gpu);

    unsigned int f_save = 100;

    OutputData output;
    output = OutputData::create("output.h5");
    output.write_scalar("density", rho);
    output.write_scalar("vel_x", vel_x);
    output.write_scalar("vel_y", vel_y);

    for (int t = 0; t < val; t++) {
        // Collision step
        collide(f1_gpu, omega);

        // Streaming step
        stream(f1_gpu, f2_gpu);

        // Decide when to save / export the data
        if ((t + 1) % f_save == 0) {
            printf("Saving data from timestep %d\n", t + 1);
            calculate_flow_properties(f2_gpu, rho_gpu, vel_x_gpu, vel_y_gpu);
            checkError(cudaMemcpy(  rho,   rho_gpu, Scalar::size*sizeof(double), cudaMemcpyDeviceToHost));
            checkError(cudaMemcpy(vel_x, vel_x_gpu, Scalar::size*sizeof(double), cudaMemcpyDeviceToHost));
            checkError(cudaMemcpy(vel_y, vel_y_gpu, Scalar::size*sizeof(double), cudaMemcpyDeviceToHost));
            output.append_scalar("density", rho);
            output.append_scalar("vel_x", vel_x);
            output.append_scalar("vel_y", vel_y);
        }

        // Swap lattices to repeat process
        double* temp = f1_gpu;
        f1_gpu = f2_gpu;
        f2_gpu = temp;
    }

    output.close();

    delete rho;
    delete vel_x;
    delete vel_y;

    checkError(cudaFree(rho_gpu));
    checkError(cudaFree(vel_x_gpu));
    checkError(cudaFree(vel_y_gpu));

    checkError(cudaFree(f1_gpu));
    checkError(cudaFree(f2_gpu));

    cudaDeviceReset();

    return 0;
}