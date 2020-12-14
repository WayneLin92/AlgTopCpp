#ifndef MYCUDA_H
#define MYCUDA_H

#define BENCHMARK_CUDA

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#ifdef BENCHMARK_CUDA
#include <iostream>
#include <chrono>
#endif

#include <vector>

constexpr int MAX_THREADS = 1024;
constexpr int MAX_INT = 0x7fffffff;

namespace cuda {
    /* Wrapper of Cuda C functions */
    void SetDevice(int device); /* Choose which GPU to run on, change this on a multi-GPU system. */
    void Malloc(void** devPtr, size_t bytes);
    void FillZero(void* devPtr, size_t bytes);
    void Memcpy(void* dst, const void* src, size_t bytes, cudaMemcpyKind kind);
    void CheckLastError();
    void DeviceSynchronize(); /* cudaDeviceSynchronize waits for the kernel to finish. */

    /* A wrapper of pointer to an array in GPU */
    class ArrayInt
    {
    public:
        ArrayInt() = default;
        explicit ArrayInt(size_t size) { Malloc((void**)&dev, size * sizeof(int)); }
        ArrayInt(const int* a, size_t size) : ArrayInt(size) { Memcpy(dev, a, size * sizeof(int), cudaMemcpyHostToDevice); }
        explicit ArrayInt(const std::vector<int>& a) : ArrayInt(a.data(), a.size()) { }
        ~ArrayInt() { cudaFree(dev); }

        void init(const std::vector<int>& a) { /* Can initialize only when dev == nullptr */
            if (dev) throw "e07509d2";
            Malloc((void**)&dev, a.size() * sizeof(int));
            Memcpy(dev, a.data(), a.size() * sizeof(int), cudaMemcpyHostToDevice);
        }
    public:
        int* data() const { return dev; }
    private:
        int* dev = nullptr;
    };

    inline void Memcpy(const ArrayInt& dst, const std::vector<int>& src) {
        Memcpy(dst.data(), src.data(), src.size() * sizeof(int), cudaMemcpyHostToDevice);
    }
    inline void Memcpy(const ArrayInt& dst, const ArrayInt& src, size_t size) {
        Memcpy(dst.data(), src.data(), size * sizeof(int), cudaMemcpyDeviceToDevice);
    }
    inline void Memcpy(std::vector<int>& dst, const ArrayInt& src, size_t size) {
        Memcpy(dst.data(), src.data(), size * sizeof(int), cudaMemcpyDeviceToHost);
    }

    /* Utilities */
    unsigned int nextPow2(unsigned int x);
    void getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads, int& blocks, int& threads);

    /* Functions using Cuda */
    void Decompress(const int* dev_in, int* dev_out, int size_in);
    inline void Decompress(const cuda::ArrayInt& dev_in, const cuda::ArrayInt& dev_out, int size_in) { Decompress(dev_in.data(), dev_out.data(), size_in); }
    
    void Sum(const int* dev_in, int* dev_out, size_t size_in);
    inline void Sum(const cuda::ArrayInt& dev_in, const cuda::ArrayInt& dev_out, size_t size_in) { Sum(dev_in.data(), dev_out.data(), size_in); }
    int Sum(const int* dev_in, size_t size_in);
    
    void Add(const int* dev_a, const int* dev_b, int* dev_out, int size);
    inline void Add(const cuda::ArrayInt& dev_a, const cuda::ArrayInt& dev_b, cuda::ArrayInt& dev_out, int size) { Add(dev_a.data(), dev_b.data(), dev_out.data(), size); }
    
    void Min(const int* dev_in, int* dev_out, size_t size_in);
    inline void Min(const cuda::ArrayInt& dev_in, const cuda::ArrayInt& dev_out, size_t size_in) { Min(dev_in.data(), dev_out.data(), size_in); }
    
    void MinIndex(const int* dev_in, int* dev_out, size_t size_in); /* Return the index of the first non-zero element. */
    inline void MinIndex(const cuda::ArrayInt& dev_in, const cuda::ArrayInt& dev_out, size_t size_in) { MinIndex(dev_in.data(), dev_out.data(), size_in); }
    int MinIndex(const int* dev_in, size_t size_in);
}

#ifdef BENCHMARK_CUDA
class Timer {
public:
    Timer() : elapsed_recent(0), bPrinted(false) { saved_time = std::chrono::system_clock::now(); }
    ~Timer() {
        if (!bPrinted) {
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed = end - saved_time;
            std::cout << "\033[0;32m" << "Elapsed time: " << elapsed.count() << "s\033[0m\n";
        }
    }
    void print(const std::string& msg = "") {
        bPrinted = true;
        auto end = std::chrono::system_clock::now();
        elapsed_recent = end - saved_time;
        std::cout << "\033[0;32m" << msg << "Elapsed time: " << elapsed_recent.count() << "s\033[0m\n";
        saved_time = std::chrono::system_clock::now();
    }
    void start() { saved_time = std::chrono::system_clock::now(); }
    std::chrono::duration<double> get_elapsed_recent() const { return elapsed_recent; }
private:
    std::chrono::time_point<std::chrono::system_clock> saved_time;
    std::chrono::duration<double> elapsed_recent;
    bool bPrinted;
};
#endif

#endif /* MYCUDA_H */
