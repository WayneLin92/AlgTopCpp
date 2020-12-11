/* This files implement cuda kernal functions */
#include "mycuda.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include <iostream>
#include <vector>

constexpr int MAX_THREADS = 1024;
#define LAUNCH_KERNEL(k, s, ...) if (s < MAX_THREADS) k<<<1, s>>>(__VA_ARGS__); else k<<<(s + MAX_THREADS - 1) / MAX_THREADS, MAX_THREADS>>>(__VA_ARGS__);

using array = std::vector<int>;
using array2d = std::vector<array>;

__global__ void InitKernel(int* a, int n)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < n)
        a[i] = 0;
}

__global__ void AddKernel(const int* a, const int* b, int* c, int n) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < n)
        c[i] = a[i] + b[i];
}

__global__ void DecompressKernel(const int* a, int* b, int offset, int size_a)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < size_a)
        b[offset + a[i]] = 1;
}

/* Return the row echelon form
** the input and out matrices are both in the form of compressed sparse rows.
*/
auto echelonCuda(const array2d& matrix_csr)
{
    cuda::SetDevice(0);
    cuda::CheckLastError();

    /* Create the sparse matrix dev_m in GPU. */
    cuda::ArrayCuda dev_m;
    std::vector<cuda::ArrayCuda> dev_rows;
    int nRows = (int)matrix_csr.size();
    int nColumns = 0;
    dev_rows.resize(nRows);
    for (int i = 0; i < nRows; ++i) {
        dev_rows[i].init(matrix_csr[i]);
        if (!matrix_csr[i].empty() && nColumns < matrix_csr[i].back() + 1)
            nColumns = matrix_csr[i].back() + 1;
    }
    int size_m = nRows * nColumns;
    dev_m.init(size_m);
    LAUNCH_KERNEL(InitKernel, size_m, dev_m.data(), size_m)
    cuda::CheckLastError();
    cuda::DeviceSynchronize();
    for (int i = 0; i < (int)matrix_csr.size(); ++i) {
        int size_rowi = (int)matrix_csr[i].size();
        LAUNCH_KERNEL(DecompressKernel, size_rowi, dev_rows[i].data(), dev_m.data(), i * nColumns, size_rowi)
        cuda::CheckLastError();
    }
    cuda::DeviceSynchronize();

    /* Reduce the rows */

    array result;
    result.resize(size_m);
    cuda::Memcpy(result.data(), dev_m.data(), size_m * sizeof(int), cudaMemcpyDeviceToHost);

    return result;
}

void addCPU(int* c, const int* a, const int* b, unsigned int size)
{
    for (size_t i = 0; i < size; ++i)
        c[i] = a[i] + b[i];
}

int runTest();

int main()
{
    return runTest();

    array2d m = { {1, 3, 4}, {1, 2, 4}, {3}, {0, 1, 2, 3, 4} };
    Timer timer;
    array sparse_m = echelonCuda(m);
    timer.print();

    /* cudaDeviceReset must be called before exiting in order for profiling and
    ** tracing tools such as Nsight and Visual Profiler to show complete traces. */
    if (cudaDeviceReset() != cudaSuccess) {
        std::cerr << "cudaDeviceReset failed!";
        return 1;
    }

    return 0;
}


