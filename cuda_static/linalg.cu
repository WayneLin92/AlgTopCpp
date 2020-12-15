#include "mycuda.h"
#include "mycuda_public.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>


namespace cuda {

    /********** Reduce the matrix `m` over F_2 by the `i`th row **********/
    __global__ void EchelonKernel(int* m, int nRows, int nColumns, int i, int j) {
        extern __shared__ int sdata[];
        int tid = threadIdx.x;
        int i1 = i + tid + 1;
        int j1 = blockIdx.x;
        if (tid == 0)
            sdata[0] = m[i * nColumns + j];
        int m_i1j = m[i1 * nColumns + j];
        __syncthreads();

        while (i1 < nRows) {
            m[i1 * nColumns + j1] ^= m[i * nColumns + j1] & sdata[0] & m_i1j;
            i1 += blockDim.x;
        }
    }

    /* Wrapper for the lauch of the kernel */
    void Echelon(cuda::ArrayInt& dev_m, int nRows, int nColumns, int i, int j)
    {
        int threads = nRows - i - 1;
        if (!threads)
            return;
        int blocks = nColumns;
        if (threads > MAX_THREADS)
            threads = MAX_THREADS;
        dim3 dimBlock(threads, 1, 1);
        dim3 dimGrid(blocks, 1, 1);
        EchelonKernel <<<dimGrid, dimBlock, sizeof(int)>>> (dev_m.data(), nRows, nColumns, i, j);
    }

    /********** replace the entries of the matrix `m` with indices for compression **********/
    __global__ void ReplaceWithIndicesKernel(int* m, int index, int size) {
        int i = threadIdx.x + blockIdx.x * blockDim.x;
        if (i < size)
            m[i] = m[i] ? index + i : -1;
    }

    /* Wrapper for the lauch of the kernel */
    void ReplaceWithIndices(int* m, int index, int size)
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
        ReplaceWithIndicesKernel<<<numBlocks, numThreads>>>(m, index, size);
#ifdef _DEBUG
        cuda::CheckLastError();
        cuda::DeviceSynchronize();
#endif
    }

    struct is_nonnegative
    {
        __host__ __device__ bool operator()(const int x) {
            return x >= 0;
        }
    };

    /* Return the row echelon form
    ** the input and out matrices are both in the form of compressed sparse rows.
    */
    array2d EchelonCuda(const array2d& matrix_csr)
    {
        /* Create the sparse matrix dev_m in GPU. */
        int nRows = (int)matrix_csr.size();
        int nColumns = 0;
        std::vector<cuda::ArrayInt> dev_rows(nRows);
        for (int i = 0; i < nRows; ++i) {
            dev_rows[i].init(matrix_csr[i]);
            if (!matrix_csr[i].empty() && nColumns < matrix_csr[i].back() + 1)
                nColumns = matrix_csr[i].back() + 1;
        }
        if (nRows == 0 || nColumns == 0)
            return array2d{};
        int size_m = nRows * nColumns;
        cuda::ArrayInt dev_m(size_m);
        cuda::FillZero(dev_m.data(), size_m * sizeof(int));
        for (int i = 0; i < (int)matrix_csr.size(); ++i) {
            int size_rowi = (int)matrix_csr[i].size();
            cuda::Decompress(dev_rows[i].data(), dev_m.data() + i * nColumns, size_rowi);
#ifdef _DEBUG
            cuda::CheckLastError();
            cuda::DeviceSynchronize();
#endif
        }

        /* Reduce the rows */
        array2d result;
        cuda::ArrayInt tmp(nColumns);
        for (int i = 0; i < nRows; ++i) {
            int index = cuda::MinIndex(dev_m.data() + i * nColumns, nColumns);
            if (index != INT_MAX) {
                Echelon(dev_m, nRows, nColumns, i, index);
#ifdef _DEBUG
                cuda::CheckLastError();
                cuda::DeviceSynchronize();
#endif
                int size_rowi = cuda::Sum(dev_m.data() + i * nColumns, nColumns);
                result.emplace_back(size_rowi);
                ReplaceWithIndices(dev_m.data() + i * nColumns + index, index, nColumns - index);
#ifdef _DEBUG
                cuda::CheckLastError();
                cuda::DeviceSynchronize();
#endif
                thrust::device_ptr<int> thrust_m = thrust::device_pointer_cast(dev_m.data());
                thrust::copy_if(thrust_m + i * nColumns, thrust_m + (i + 1) * nColumns, tmp.data(), is_nonnegative());
                cuda::Memcpy(result.back(), tmp, result.back().size());
            }
        }
        
        return result;
    }
}