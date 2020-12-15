#include "mycuda.h"
#include "mycuda_public.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

constexpr int WARP_SIZE = 32;

/* perform first level of reduction */
template <unsigned int blockSize>
__global__ void ReduceAddKernel(const int* __restrict__ g_idata, int* __restrict__ g_odata, unsigned int n)
{
    extern __shared__ int sdata[];
    unsigned int tid = threadIdx.x;
    unsigned int gridSize2 = blockSize * gridDim.x * 2;

    /* we reduce multiple elements per thread.  The number is determined by gridSize. 
    ** More blocks will result in a larger gridSize and therefore fewer elements per thread */
    unsigned int i = blockIdx.x * blockSize * 2 + threadIdx.x;
    int mySum = 0;
    while (i < n) {
        mySum += g_idata[i];
        if ((i + blockSize) < n)
            mySum += g_idata[i + blockSize];
        i += gridSize2;
    }

    /* Reduce within warp using __shfl_down_sync */
    for (int offset = min(blockSize, WARP_SIZE) / 2; offset > 0; offset /= 2)
        mySum += __shfl_down_sync(0xffffffff, mySum, offset);
    if ((tid % WARP_SIZE) == 0) /* each warp puts its local sum into shared memory */
        sdata[tid / WARP_SIZE] = mySum;

    __syncthreads();

    /* Reduce shared memory using __shfl_down_sync  */
    const unsigned int size_share_memory = (blockSize / WARP_SIZE) > 0 ? (blockSize / WARP_SIZE) : 1; /* size_share_memory <= 1024/32=32 */
    const unsigned int mask_ballot = __ballot_sync(0xffffffff, tid < size_share_memory);
    if (tid < size_share_memory) {
        mySum = sdata[tid];
        for (int offset = size_share_memory / 2; offset > 0; offset /= 2)
            mySum += __shfl_down_sync(mask_ballot, mySum, offset);
    }

    /* write result for this block to global mem */
    if (tid == 0)
        g_odata[blockIdx.x] = mySum;
}

/* Reduce the array `dev_in` of size `size` to the array `dev_out` of size `blocks` */
void ReduceAdd(int size, int threads, int blocks, const int* dev_in, int* dev_out)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    int smemSize = ((threads / WARP_SIZE) + 1) * sizeof(int);
    switch (threads) {
    case 1024:
        ReduceAddKernel<1024><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 512:
        ReduceAddKernel<512> <<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 256:
        ReduceAddKernel<256> <<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 128:
        ReduceAddKernel<128> <<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 64:
        ReduceAddKernel<64>  <<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 32:
        ReduceAddKernel<32>  <<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 16:
        ReduceAddKernel<16>  <<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case  8:
        ReduceAddKernel<8>   <<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case  4:
        ReduceAddKernel<4>   <<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case  2:
        ReduceAddKernel<2>   <<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case  1:
        ReduceAddKernel<1>   <<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    }
}

unsigned int cuda::nextPow2(unsigned int x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

/* Compute numbers of blocks and threads for reduction algorithms */
void cuda::getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads, int& blocks, int& threads) {
    threads = (n < maxThreads * 2) ? nextPow2((n + 1) / 2) : maxThreads;
    blocks = (n + (threads * 2 - 1)) / (threads * 2);
    if (blocks > maxBlocks)
        blocks = maxBlocks;
}

void cuda::Sum(const int* dev_in, int* dev_out, size_t size_in) /* size_in should be positive */
{
    int maxThreads = 256;
    int maxBlocks = 64;
    int numBlocks = 0;
    int numThreads = 0;
    getNumBlocksAndThreads((int)size_in, maxBlocks, maxThreads, numBlocks, numThreads);
    cuda::ArrayInt dev_c(numBlocks), dev_tmp(numBlocks);
    ReduceAdd((int)size_in, numThreads, numBlocks, dev_in, dev_c.data());
#ifdef _DEBUG
    CheckLastError();
    DeviceSynchronize();
#endif

    int s = numBlocks;
    while (s > 1) {
        int threads = 0, blocks = 0;
        getNumBlocksAndThreads(s, maxBlocks, maxThreads, blocks, threads);
        Memcpy(dev_tmp, dev_c, s);
        ReduceAdd(s, threads, blocks, dev_tmp.data(), dev_c.data());
#ifdef _DEBUG
        CheckLastError();
        DeviceSynchronize();
#endif
        s = (s + (threads * 2 - 1)) / (threads * 2);
    }
    Memcpy(dev_out, dev_c.data(), sizeof(int), cudaMemcpyDeviceToDevice);
}

int cuda::Sum(const int* dev_in, size_t size_in)
{
    cuda::ArrayInt dev_out(1);
    Sum(dev_in, dev_out.data(), size_in);
    int sum;
     cuda::Memcpy(&sum, dev_out.data(), sizeof(int), cudaMemcpyDeviceToHost);
    return sum;
}