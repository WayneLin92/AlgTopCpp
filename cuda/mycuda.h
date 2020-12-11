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

namespace cuda {
    void SetDevice(int device); /* Choose which GPU to run on, change this on a multi-GPU system. */
    void Malloc(void** devPtr, size_t size);
    void Memcpy(void* dst, const void* src, size_t count, cudaMemcpyKind kind);
    void CheckLastError();
    void DeviceSynchronize(); /* cudaDeviceSynchronize waits for the kernel to finish. */

    /* A wrapper of pointer to an array in GPU */
    class ArrayCuda
    {
    public:
        ArrayCuda() : dev(nullptr) {}
        ~ArrayCuda() { cudaFree(dev); }
        void init(size_t n) { Malloc((void**)&dev, n * sizeof(int)); }
        void init(const int* a, size_t n) { init(n); Memcpy(dev, a, n * sizeof(int), cudaMemcpyHostToDevice); }
        void init(const std::vector<int>& a) { init(a.data(), a.size()); }
    public:
        int* data() const { return dev; }
    private:
        int* dev;
    };
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
    std::chrono::duration<double> get_elapsed_recent() const { return elapsed_recent; }
private:
    std::chrono::time_point<std::chrono::system_clock> saved_time;
    std::chrono::duration<double> elapsed_recent;
    bool bPrinted;
};
#endif

#endif
