#include "mycuda.h"
#include "mycuda_public.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

constexpr int WARP_SIZE = 32;

/* perform first level of reduction */
template <unsigned int blockSize>
__global__ void ReduceMinIndexKernel(const int* __restrict__ g_idata, int* __restrict__ g_odata, unsigned int n)
{
    extern __shared__ int sdata[];
    unsigned int tid = threadIdx.x;
    unsigned int gridSize2 = blockSize * gridDim.x * 2;

    /* we reduce multiple elements per thread.  The number is determined by gridSize.
    ** More blocks will result in a larger gridSize and therefore fewer elements per thread */
    unsigned int i = blockIdx.x * blockSize * 2 + threadIdx.x;
    int myIndex = INT_MAX;
    while (i < n) {
        myIndex = min(myIndex, g_idata[i] > 0 ? i : INT_MAX);
        if ((i + blockSize) < n)
            myIndex = min(myIndex, g_idata[i + blockSize] > 0 ? i + blockSize : INT_MAX);
        i += gridSize2;
    }

    /* Reduce within warp using __shfl_down_sync */
    for (int offset = min(blockSize, WARP_SIZE) / 2; offset > 0; offset /= 2)
        myIndex = min(myIndex, __shfl_down_sync(0xffffffff, myIndex, offset));
    if ((tid % WARP_SIZE) == 0) /* each warp puts its local sum into shared memory */
        sdata[tid / WARP_SIZE] = myIndex;
    __syncthreads();

    /* Reduce shared memory using __shfl_down_sync  */
    const unsigned int size_share_memory = (blockSize / WARP_SIZE) > 0 ? (blockSize / WARP_SIZE) : 1; /* size_share_memory <= 1024/32=32 */
    const unsigned int mask_ballot = __ballot_sync(0xffffffff, tid < size_share_memory);
    if (tid < size_share_memory) {
        myIndex = sdata[tid];
        for (int offset = size_share_memory / 2; offset > 0; offset /= 2)
            myIndex = min(myIndex, __shfl_down_sync(mask_ballot, myIndex, offset));
    }

    /* write result for this block to global mem */
    if (tid == 0)
        g_odata[blockIdx.x] = myIndex;
}

/* Wrapper for the lauch of the kernel
** Reduce the array `dev_in` of size `size` to the array `dev_out` of size `blocks` */
void ReduceMinIndex(const int* dev_in, int* dev_out, int threads, int blocks, int size)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    int smemSize = ((threads / WARP_SIZE) + 1) * sizeof(int);
    switch (threads) {
    case 1024:
        ReduceMinIndexKernel<1024><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 512:
        ReduceMinIndexKernel<512><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 256:
        ReduceMinIndexKernel<256><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 128:
        ReduceMinIndexKernel<128><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 64:
        ReduceMinIndexKernel<64><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 32:
        ReduceMinIndexKernel<32><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case 16:
        ReduceMinIndexKernel<16><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case  8:
        ReduceMinIndexKernel<8><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case  4:
        ReduceMinIndexKernel<4><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case  2:
        ReduceMinIndexKernel<2><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    case  1:
        ReduceMinIndexKernel<1><<<dimGrid, dimBlock, smemSize>>>(dev_in, dev_out, size);
        break;
    }
}

void ReduceMin(int size, int threads, int blocks, const int* dev_in, int* dev_out);

void cuda::MinIndex(const int* dev_in, int* dev_out, size_t size_in)
{
    if (size_in == 0) {
        *dev_out = INT_MAX;
        return;
    }
    int maxThreads = 256;
    int maxBlocks = 64;
    int numBlocks = 0;
    int numThreads = 0;
    getNumBlocksAndThreads((int)size_in, maxBlocks, maxThreads, numBlocks, numThreads);
    ArrayInt dev_c(numBlocks), dev_tmp(numBlocks);
    ReduceMinIndex(dev_in, dev_c.data(), numThreads, numBlocks, (int)size_in);
#ifdef _DEBUG
    CheckLastError();
    DeviceSynchronize();
#endif

    int s = numBlocks;
    while (s > 1) {
        int threads = 0, blocks = 0;
        getNumBlocksAndThreads(s, maxBlocks, maxThreads, blocks, threads);
        Memcpy(dev_tmp, dev_c, s);
        ReduceMin(s, threads, blocks, dev_tmp.data(), dev_c.data());
#ifdef _DEBUG
        CheckLastError();
        DeviceSynchronize();
#endif
        s = (s + (threads * 2 - 1)) / (threads * 2);
    }
    Memcpy(dev_out, dev_c.data(), sizeof(int), cudaMemcpyDeviceToDevice);
}

int cuda::MinIndex(const int* dev_in, size_t size_in)
{
    cuda::ArrayInt dev_out(1);
    MinIndex(dev_in, dev_out.data(), size_in);
    int index;
    cuda::Memcpy(&index, dev_out.data(), sizeof(int), cudaMemcpyDeviceToHost);
    return index;
}
