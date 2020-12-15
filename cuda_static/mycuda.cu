#include "mycuda.h"
#include <iostream>

void cuda::SetDevice(int device) {
    if (cudaSetDevice(device) != cudaSuccess) {
        std::cerr << "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?";
        throw "e31e20a";
    }
}

void cuda::Malloc(void** devPtr, size_t bytes)
{
    cudaError_t cudaStatus = cudaMalloc(devPtr, bytes);
    if (cudaStatus != cudaSuccess) {
        std::cerr << "cudaMalloc failed: " << cudaGetErrorString(cudaStatus) << '\n';
        throw "1e32d383";
    }
}

void cuda::FillZero(void* devPtr, size_t bytes)
{
    cudaError_t cudaStatus = cudaMemset(devPtr, 0, bytes);
    if (cudaStatus != cudaSuccess) {
        std::cerr << "FillZero failed!: " << cudaGetErrorString(cudaStatus) << '\n';
        throw "798e95c";
    }
}

void cuda::Memcpy(void* dst, const void* src, size_t bytes, cudaMemcpyKind kind)
{
    cudaError_t cudaStatus = cudaMemcpy(dst, src, bytes, kind);
    if (cudaStatus != cudaSuccess) {
        std::cerr << "cudaMemcpy failed: " << cudaGetErrorString(cudaStatus) << '\n';
        throw "ac621fd2";
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
        std::cerr << "cudaDeviceSynchronize failed: " << cudaGetErrorString(cudaStatus) << '\n';
        throw "19149f41";
    }
}