#include "database.h"
#include <chrono>

int main_test(int argc, char** argv)
{
	array2d spaceV = { {1, 2, 3, 4}, {2, 3, 4}, {3, 4}, {4} };
	array2d spaceW = { {2, 3}, {5} };
	bool bException = false;
	try {
		array2d quotient = quotient_space(spaceV, spaceW);
		std::cout << quotient << '\n';
	}
	catch (const char* e) {
		std::cout << e << '\n';
		bException = true;
	}
	return 0;
}