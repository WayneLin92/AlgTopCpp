#ifndef MAIN_H
#define MAIN_H

#include "database.h"
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

void generate_basis(const Database& db, const std::string& table_prefix, int t_max, bool drop_existing = false);
void generate_mon_diffs(const Database& db, const std::string& table_prefix, int r);
void generate_ss(const Database& db, const std::string& table_name_basis, const std::string& table_ss, int r);
void generate_next_page(const Database& db, const std::string& table_prefix, const std::string& table_H_prefix, int r);
int main_test(int argc, char** argv);


#endif
