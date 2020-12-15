#include "mycuda.h"
#include "mycuda_public.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

using array = std::vector<int>;
using array2d = std::vector<array>;

__global__ void TmpKernel(const int* in, int* out, int size) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < size)
        printf("i=%d, in[i]=%d\n", i, in[i]);
}

void Tmp(const cuda::ArrayInt& dev_in, cuda::ArrayInt& dev_out, int size)
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
    TmpKernel <<< numBlocks, numThreads >>> (dev_in.data() + 3, dev_out.data(), size - 3);
#ifdef _DEBUG
    cuda::CheckLastError();
    cuda::DeviceSynchronize();
#endif
}

int testTmp()
{
    int size = 24;
    std::vector<int> a(size);
    for (int i = 0; i < size; ++i)
        a[i] = i;
    std::vector<int> b(size);

    cuda::ArrayInt dev_a(a), dev_b(size);
    Tmp(dev_a, dev_b, size);
    cuda::Memcpy(b, dev_b, size);

    /*for (int i = 0; i < size; ++i)
        std::cout << "b" << i << '=' << b[i] << '\n';*/
    return 0;
}

int test_sum()
{
    array b(1);

    array a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    cuda::ArrayInt dev_a(a), dev_b(1);
    cuda::Sum(dev_a, dev_b, a.size());
    cuda::Memcpy(b, dev_b, 1);
    std::cout << b[0] << '\n';

    array a2(1030, 1);
    cuda::ArrayInt dev_a2(a2);
    cuda::Sum(dev_a2, dev_b, a2.size());
    cuda::Memcpy(b, dev_b, 1);
    std::cout << b[0] << '\n';

    array a3 = { 6, 5, 4, 5, 7, 8, 10 };
    cuda::ArrayInt dev_a3(a3);
    int sum3 = cuda::Sum(dev_a3.data() + 1, a3.size() - 1);
    std::cout << sum3 << '\n'; /* 39 */

    return 0;
}

int test_min()
{
    array b(1);
    cuda::ArrayInt dev_b(1);

    array a(128);
    for (int i = 0; i < (int)a.size(); ++i)
        a[i] = int(a.size()) - i + 1000000;
    cuda::ArrayInt dev_a(a);
    cuda::Min(dev_a, dev_b, a.size());
    cuda::Memcpy(b, dev_b, 1);
    std::cout << b[0] << '\n';

    array a2(64, 9);
    a2[31] = 7;
    a2[33] = 10;
    cuda::ArrayInt dev_a2(a2);
    cuda::Min(dev_a2, dev_b, a2.size());
    cuda::Memcpy(b, dev_b, 1);
    std::cout << b[0] << '\n';
    cuda::CheckLastError();

    array a3 = { 6, 5, 4, 5, 13, 7, 8, 10, 2 };
    cuda::ArrayInt dev_a3(a3);
    cuda::Min(dev_a3, dev_b, a3.size());
    cuda::Memcpy(b, dev_b, 1);
    std::cout << b[0] << '\n';
    cuda::CheckLastError();

    return 0;
}

int test_min_index()
{
    array b(1);
    cuda::ArrayInt dev_b(1);

    array a(128, 0);
    a[50] = 1;
    cuda::ArrayInt dev_a(a);
    cuda::MinIndex(dev_a, dev_b, a.size());
    cuda::Memcpy(b, dev_b, 1);
    std::cout << b[0] << " Expected: " << 50 << '\n';

    array a2(60, 0);
    a2[31] = 7;
    a2[33] = 10;
    cuda::ArrayInt dev_a2(a2);
    cuda::MinIndex(dev_a2, dev_b, a2.size());
    cuda::Memcpy(b, dev_b, 1);
    std::cout << b[0] << " Expected: " << 31 << '\n';
    cuda::CheckLastError();

    array a3 = { 0, 0, 1, 2, 3, 4, 5 };
    cuda::ArrayInt dev_a3(a3);
    cuda::MinIndex(dev_a3, dev_b, a3.size());
    cuda::Memcpy(b, dev_b, 1);
    std::cout << b[0] << " Expected: " << 2 << '\n';
    cuda::CheckLastError();

    return 0;
}

int test_echelon()
{
    Timer timer;

    array2d m = { {1, 3, 4}, {1, 2, 4}, {3}, {0, 1, 2, 3, 4} };
    array2d rref = cuda::EchelonCuda(m);
    for (size_t i = 0; i < rref.size(); ++i) {
        for (size_t j = 0; j < rref[i].size(); ++j)
            std::cout << rref[i][j] << ", ";
        std::cout << '\n';
    }
    std::cout << '\n';
    return 0;
}

int main()
{
    int return_code = test_echelon();

    /* cudaDeviceReset must be called before exiting in order for profiling and
    ** tracing tools such as Nsight and Visual Profiler to show complete traces. */
    if (cudaDeviceReset() != cudaSuccess) {
        std::cerr << "cudaDeviceReset failed!";
        return 1;
    }

    return return_code;
}