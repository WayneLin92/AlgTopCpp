#include "mycuda.h"
#include "mycuda_public.h"
#include "myio.h"
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
    if (size < cuda::MAX_THREADS) {
        numBlocks = 1;
        numThreads = size;
    }
    else {
        numBlocks = (size + cuda::MAX_THREADS - 1) / cuda::MAX_THREADS;
        numThreads = cuda::MAX_THREADS;
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

int test_echelon() /* Test with random matrix */
{
    //Timer timer;

    array2d m = { {84}, {143}, {84}, {143}, {137}, {158}, {143}, {63}, {14, 52, 82}, {60, 61, 82}, {49, 52, 60}, {137}, {116}, {54}, {55}, {156}, {142}, {22}, {80}, {}, {59, 62, 70}, {123}, {142}, {80, 110, 112, 141}, {46, 81, 142}, {25}, {}, {14, 49, 61}, {13, 15, 62, 113}, {54}, {55, 116}, {29}, {123}, {48, 123}, {116, 123}, {54}, {55}, {29}, {20, 30, 109}, {29}, {23, 30, 119}, {109, 118, 119}, {42, 44, 107}, {19, 43, 45}, {44, 57}, {45, 58}, {107}, {108, 117, 129}, {17}, {20, 23, 118}, {42, 57}, {43, 58}, {136}, {19}, {}, {}, {173, 182}, {128, 149}, {184, 185}, {147, 149}, {151, 155}, {34, 127}, {}, {}, {89, 96}, {177, 183}, {}, {139, 160}, {79, 106, 155}, {68, 69, 79}, {125, 127}, {128, 147}, {171, 183}, {170, 184}, {}, {159, 160}, {}, {72, 96}, {94, 96}, {}, {}, {}, {}, {41, 106, 151}, {77, 78}, {34, 125}, {87, 97}, {171, 177}, {173, 177, 180}, {145, 150}, {}, {}, {41}, {}, {148, 167}, {}, {}, {}, {}, {163, 167}, {165, 167}, {124, 140}, {72, 89}, {89, 94}, {28, 122}, {134, 135}, {139, 159}, {169, 175, 176, 179}, {40}, {}, {}, {65, 88}, {162, 168}, {101, 122}, {}, {}, {86}, {148, 163}, {148, 165}, {}, {}, {71, 93}, {72, 94}, {27, 31}, {168, 174}, {161, 172}, {}, {163, 165}, {38, 39, 100, 104}, {}, {}, {}, {}, {86}, {86}, {}, {28, 37, 76, 101}, {16, 53, 115}, {75, 76}, {135, 153}, {26, 120}, {126, 144}, {126, 146}, {132, 162}, {162, 174}, {}, {99, 120}, {33}, {}, {144, 146}, {64, 66}, {66, 74, 85, 91}, {85, 90, 91}, {37, 75}, {}, {132, 174}, {133, 154}, {11, 64, 154}, {64, 74, 90}, {26, 35, 67, 99}, {12, 35} };
    array2d rref = cuda::EchelonCuda(m);
    
    std::cout << rref << '\n';
    return 0;
}

int main1()
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