#include "mycuda.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

using T = int;
constexpr int WARP_SIZE = 32;

template <unsigned int blockSize>
__global__ void reduceKernel(const T* __restrict__ g_idata, T* __restrict__ g_odata, unsigned int n)
{
    /* perform first level of reduction,
    ** reading from global memory, writing to shared memory */

    extern __shared__ int sdata[];
    unsigned int tid = threadIdx.x;
    unsigned int gridSize2 = blockSize * gridDim.x * 2;
    unsigned int maskLength = blockSize & (WARP_SIZE - 1);
    maskLength = (maskLength > 0) ? (WARP_SIZE - maskLength) : maskLength; // ?
    const unsigned int mask = (0xffffffff) >> maskLength;

    /* we reduce multiple elements per thread.  The number is determined by gridSize. 
    ** More blocks will result in a larger gridSize and therefore fewer elements per thread */
    unsigned int i = blockIdx.x * blockSize * 2 + threadIdx.x;
    T mySum = 0;
    while (i < n) {
        mySum += g_idata[i];
        if ((i + blockSize) < n)
            mySum += g_idata[i + blockSize];
        i += gridSize2;
    }

    /* Reduce within warp using __shfl_down_sync */
    for (int offset = WARP_SIZE / 2; offset > 0; offset /= 2)
        mySum += __shfl_down_sync(mask, mySum, offset);
    if ((tid % WARP_SIZE) == 0) /* each warp puts its local sum into shared memory */
        sdata[tid / WARP_SIZE] = mySum;

    __syncthreads();

    /* Reduce shared memory */
    const unsigned int shmem_extent = (blockSize / WARP_SIZE) > 0 ? (blockSize / WARP_SIZE) : 1;
    const unsigned int mask_ballot = __ballot_sync(mask, tid < shmem_extent);
    if (tid < shmem_extent) {
        mySum = sdata[tid];
        /* Reduce within warp using __shfl_down_sync */
        for (int offset = WARP_SIZE / 2; offset > 0; offset /= 2)
            mySum += __shfl_down_sync(mask_ballot, mySum, offset);
    }

    // write result for this block to global mem
    if (tid == 0)
        g_odata[blockIdx.x] = mySum;
}

void reduce(int size, int threads, int blocks, T* d_idata, T* d_odata)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    int smemSize = ((threads / WARP_SIZE) + 1) * sizeof(T);
    switch (threads) {
    case 1024:
        reduceKernel<1024><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    case 512:
        reduceKernel<512> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    case 256:
        reduceKernel<256> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    case 128:
        reduceKernel<128> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    case 64:
        reduceKernel<64>  <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    case 32:
        reduceKernel<32>  <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    case 16:
        reduceKernel<16>  <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    case  8:
        reduceKernel<8>   <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    case  4:
        reduceKernel<4>   <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    case  2:
        reduceKernel<2>   <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    case  1:
        reduceKernel<1>   <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    }
}

unsigned int nextPow2(unsigned int x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

void getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads, int& blocks, int& threads) {
    threads = (n < maxThreads * 2) ? nextPow2((n + 1) / 2) : maxThreads;
    blocks = (n + (threads * 2 - 1)) / (threads * 2);
    if (blocks > maxBlocks)
        blocks = maxBlocks;
}

int runTest()
{

    int size = 1 << 24;    /* number of elements to reduce */
    int maxThreads = 256;  /* number of threads per block */
    int maxBlocks = 64;

    /* create random input data on CPU */
    unsigned int bytes = size * sizeof(T);
    T* h_idata = (T*)malloc(bytes);
    for (int i = 0; i < size; i++)
        h_idata[i] = (T)(rand() & 0xFF);  /* Keep the numbers small so we don't get truncation error in the sum */

    Timer timer_cpu;
    int cpu_result = 0;
    for (int i = 0; i < size; i++)
        cpu_result += h_idata[i];
    timer_cpu.print("CPU: ");

    int numBlocks = 0;
    int numThreads = 0;
    getNumBlocksAndThreads(size, maxBlocks, maxThreads, numBlocks, numThreads);

    T* d_idata = NULL;
    T* d_odata = NULL;
    cuda::Malloc((void**)&d_idata, bytes);
    cuda::Malloc((void**)&d_odata, numBlocks * sizeof(T));

    cuda::Memcpy(d_idata, h_idata, bytes, cudaMemcpyHostToDevice);
    cuda::Memcpy(d_odata, h_idata, numBlocks * sizeof(T), cudaMemcpyHostToDevice);

    // warm-up
    reduce(size, numThreads, numBlocks, d_idata, d_odata);
    cuda::CheckLastError();
    cuda::DeviceSynchronize();

    Timer timer_gpu;
    int gpu_result = 0;
    for (int i = 0; i < 100; ++i) {

        T* d_intermediateSums;
        cuda::Malloc((void**)&d_intermediateSums, numBlocks * sizeof(T));

        reduce(size, numThreads, numBlocks, d_idata, d_odata);
        cuda::CheckLastError();

        int s = numBlocks;
        while (s > 1)
        {
            int threads = 0, blocks = 0;
            getNumBlocksAndThreads(s, maxBlocks, maxThreads, blocks, threads);
            cuda::Memcpy(d_intermediateSums, d_odata, s * sizeof(T), cudaMemcpyDeviceToDevice);
            reduce(s, threads, blocks, d_intermediateSums, d_odata);
            s = (s + (threads * 2 - 1)) / (threads * 2);
        }

        cuda::DeviceSynchronize();
        cuda::Memcpy(&gpu_result, d_odata, sizeof(T), cudaMemcpyDeviceToHost);
    }
    timer_gpu.print("GPU * 100: ");

    std::cout << "Result (CPU) = " << cpu_result << '\n';
    std::cout << "Result (GPU) = " << gpu_result << '\n';
    std::cout << "Time (GPU) = " << timer_gpu.get_elapsed_recent().count() << '\n';
    double reduceTime = timer_gpu.get_elapsed_recent().count() / 100;
    printf("Reduction, Throughput = %.4f GB/s, Time = %.5f s, Size = %u Elements, NumDevsUsed = %d, Workgroup = %u\n",
        1.0e-9 * ((double)bytes) / reduceTime, reduceTime, size, 1, numThreads);

    return cpu_result == gpu_result ? 0 : 1;
}