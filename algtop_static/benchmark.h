#pragma once

#include <chrono>

class Timer {
public:
	Timer() : bPrinted(false) { saved_time = std::chrono::system_clock::now(); }
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
		std::chrono::duration<double> elapsed = end - saved_time;
		std::cout << "\033[0;32m" << msg << "Elapsed time: " << elapsed.count() << "s\033[0m\n";
		saved_time = end;
	}
private:
	std::chrono::time_point<std::chrono::system_clock> saved_time;
	bool bPrinted;
};