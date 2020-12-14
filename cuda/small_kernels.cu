/* This files implement cuda kernal functions */
#include "mycuda.h"
#include "mycuda_public.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include <iostream>
#include <vector>


/********** Add two vectors over F_2 **********/
__global__ void AddKernel(const int* a, const int* b, int* out, int size) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < size)
        out[i] = a[i] ^ b[i];
}

void cuda::Add(const int* dev_a, const int* dev_b, int* dev_out, int size)
{
    int numBlocks;
    int numThreads;
    if (size < MAX_THREADS) {
        numBlocks = 1;
        numThreads = size;
    }
    else {
        numBlocks = (size + MAX_THREADS - 1) / MAX_THREADS;
        numThreads = MAX_THREADS;
    }
    AddKernel<<<numBlocks, numThreads>>>(dev_a, dev_b, dev_out, size);
#ifdef _DEBUG
    cuda::CheckLastError();
    cuda::DeviceSynchronize();
#endif
}

/********** Decompress a vector **********/
__global__ void DecompressKernel(const int* in, int* out, int size_in)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < size_in)
        out[in[i]] = 1;
}

void cuda::Decompress(const int* dev_in, int* dev_out, int size_in)
{
    int numBlocks;
    int numThreads;
    if (size_in < MAX_THREADS) {
        numBlocks = 1;
        numThreads = size_in;
    }
    else {
        numBlocks = (size_in + MAX_THREADS - 1) / MAX_THREADS;
        numThreads = MAX_THREADS;
    }
    DecompressKernel<<<numBlocks, numThreads>>>(dev_in, dev_out, size_in);
#ifdef _DEBUG
    cuda::CheckLastError();
    cuda::DeviceSynchronize();
#endif
}


