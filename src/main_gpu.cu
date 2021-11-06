#include <array>
#include <chrono>
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

void DisplayHeader()
{
    const int kb = 1024;
    const int mb = kb * kb;
    cout << "\nCUDA version:   v" << CUDART_VERSION << endl;

    int devCount;
    cudaGetDeviceCount(&devCount);
    cout << "CUDA Devices: " << endl << endl;

    for(int i = 0; i < devCount; ++i)
    {
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);
        cout << i << ": " << props.name << ": " << props.major << "." << props.minor << endl;
        cout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << endl;
        cout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << endl;
        cout << "  Constant memory: " << props.totalConstMem / kb << "kb" << endl;
        cout << "  Block registers: " << props.regsPerBlock << endl << endl;

        cout << "  Warp size:         " << props.warpSize << endl;
        cout << "  Threads per block: " << props.maxThreadsPerBlock << endl;
        cout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1]  << ", " << props.maxThreadsDim[2] << " ]" << endl;
        cout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", " << props.maxGridSize[1]  << ", " << props.maxGridSize[2] << " ]" << endl;
        cout << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Expected filename" << endl;
        return -1;
    }

    DisplayHeader();

    double nu_lb = 0.092;                     // Lattice dynamic viscosity
    double omega = 1.0 / (3. * nu_lb + 0.5);  // Relaxation parameter

    InputData config = InputData::open(argv[1]);

    Field::set_field_shape(config.shape_x, config.shape_y, Q);
    Scalar::set_scalar_shape(config.shape_x, config.shape_y);

    size_x = config.shape_x;
    size_y = config.shape_y;

    set_gpu_size<<<1, 1>>>(size_x, size_y);

    int* type_lattice = new int[Scalar::size];

    config.configure_domain(type_lattice);

    double* rho   = new double[Scalar::size];
    double* vel_x = new double[Scalar::size];
    double* vel_y = new double[Scalar::size];

    double* rho_gpu;
    double* vel_x_gpu;
    double* vel_y_gpu;

    cudaEvent_t begin, end;
    checkError(cudaEventCreate(&begin));
    checkError(cudaEventCreate(&end));

    checkError(cudaMalloc((void**)&rho_gpu,   Scalar::size*sizeof(double)));
    checkError(cudaMalloc((void**)&vel_x_gpu, Scalar::size*sizeof(double)));
    checkError(cudaMalloc((void**)&vel_y_gpu, Scalar::size*sizeof(double)));

    config.read_scalar("initial/density", rho);
    config.read_scalar("initial/vel_x", vel_x);
    config.read_scalar("initial/vel_y", vel_y);

    checkError(cudaMemcpy(  rho_gpu,   rho, Scalar::size*sizeof(double), cudaMemcpyHostToDevice));
    checkError(cudaMemcpy(vel_x_gpu, vel_x, Scalar::size*sizeof(double), cudaMemcpyHostToDevice));
    checkError(cudaMemcpy(vel_y_gpu, vel_y, Scalar::size*sizeof(double), cudaMemcpyHostToDevice));

    config.close();

    double* f1_gpu;
    double* f2_gpu;

    checkError(cudaMalloc((void**)&f1_gpu, Field::size*sizeof(double)));
    checkError(cudaMalloc((void**)&f2_gpu, Field::size*sizeof(double)));

    // Set up equilibrium values for initial flow
    init_eq(f1_gpu, rho_gpu, vel_x_gpu, vel_y_gpu);

    OutputData output;
    output = OutputData::create("output.h5");

    auto saveTime = chrono::microseconds::zero();
    auto startSave = chrono::high_resolution_clock::now();
    output.write_scalar("density", rho);
    output.write_scalar("vel_x", vel_x);
    output.write_scalar("vel_y", vel_y);
    auto endSave = chrono::high_resolution_clock::now();
    saveTime += chrono::duration_cast<chrono::microseconds>(endSave - startSave);

    auto start = chrono::high_resolution_clock::now();
    checkError(cudaEventRecord(begin, 0));

    for (int t = 0; t < config.timesteps; t++) {
        // Collision step
        collide(f1_gpu, omega);

        // Streaming step
        stream(f1_gpu, f2_gpu);

        // Decide when to save / export the data
        if ((config.savestep > 0) && ((t + 1) % config.savestep == 0)) {
            printf("Saving data from timestep %d\n", t + 1);
            calculate_flow_properties(f2_gpu, rho_gpu, vel_x_gpu, vel_y_gpu);
            auto startSave = chrono::high_resolution_clock::now();
            checkError(cudaMemcpy(  rho,   rho_gpu, Scalar::size*sizeof(double), cudaMemcpyDeviceToHost));
            checkError(cudaMemcpy(vel_x, vel_x_gpu, Scalar::size*sizeof(double), cudaMemcpyDeviceToHost));
            checkError(cudaMemcpy(vel_y, vel_y_gpu, Scalar::size*sizeof(double), cudaMemcpyDeviceToHost));
            output.append_scalar("density", rho);
            output.append_scalar("vel_x", vel_x);
            output.append_scalar("vel_y", vel_y);
            auto endSave = chrono::high_resolution_clock::now();
            saveTime += chrono::duration_cast<chrono::microseconds>(endSave - startSave);
        }

        // Swap lattices to repeat process
        double* temp = f1_gpu;
        f1_gpu = f2_gpu;
        f2_gpu = temp;
    }

    checkError(cudaEventRecord(end, 0));
    checkError(cudaEventSynchronize(end));
    float gpu_milliseconds = 0.0f;
    checkError(cudaEventElapsedTime(&gpu_milliseconds, begin, end));

    auto stop = chrono::high_resolution_clock::now();
    double runtime = chrono::duration_cast<chrono::milliseconds> (stop - start).count() / 1000.0;
    double gpu_runtime = 0.001*gpu_milliseconds;

    size_t nodes_updated = config.timesteps * size_t(size_x * size_y);
    size_t nodes_saved   = config.timesteps / config.savestep * size_t(size_x * size_y);

    // calculate speed in million lattice updates per second
    double speed = nodes_updated / (1E+06 * runtime);
    // calculate memory access rate in GiB/s
    double bytesPerGiB = 1024.0 * 1024.0 * 1024.0;
    double bandwidth = (
        nodes_updated * Q * 2
        + nodes_saved * 3
    ) * sizeof(double) / (runtime * bytesPerGiB);

    printf(" ----- performance information -----\n");
    printf("        timesteps: %u\n", config.timesteps);
    printf("    clock runtime: %.3f (s)\n", runtime);
    printf("      gpu runtime: %.3f (s)\n", gpu_runtime);
    printf("      output time: %.3f (ms)\n", saveTime.count() / 1000.0);
    printf("            speed: %.2f (Mlups)\n",speed);
    printf("        bandwidth: %.1f (GiB/s)\n",bandwidth);

    // Close resources

    output.close();

    delete rho;
    delete vel_x;
    delete vel_y;

    checkError(cudaEventDestroy(begin));
    checkError(cudaEventDestroy(end));

    checkError(cudaFree(rho_gpu));
    checkError(cudaFree(vel_x_gpu));
    checkError(cudaFree(vel_y_gpu));

    checkError(cudaFree(f1_gpu));
    checkError(cudaFree(f2_gpu));

    cudaDeviceReset();

    return 0;
}