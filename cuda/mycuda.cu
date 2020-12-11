#include "mycuda.h"
#include <iostream>

void cuda::SetDevice(int device) {
    if (cudaSetDevice(device) != cudaSuccess) {
        std::cerr << "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?";
        throw "e31e20a";
    }
}

void cuda::Malloc(void** devPtr, size_t size)
{
    CheckLastError();
    if (cudaMalloc(devPtr, size) != cudaSuccess) {
        std::cerr << "cudaMalloc failed!";
        throw "1e32d383";
    }
}

void cuda::Memcpy(void* dst, const void* src, size_t count, cudaMemcpyKind kind)
{
    if (cudaMemcpy(dst, src, count, kind) != cudaSuccess) {
        std::cerr << "cudaMemcpy failed!";
        throw "798e95c";
    }
}

void cuda::CheckLastError() {
    cudaError_t cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        std::cerr << "Cuda failed: " << cudaGetErrorString(cudaStatus) << '\n';
        throw "e502018d";
    }
}

void cuda::DeviceSynchronize()
{
    cudaError_t cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        std::cerr << "cudaDeviceSynchronize returned error code " << cudaStatus << '\n';
        throw "19149f41";
    }
}