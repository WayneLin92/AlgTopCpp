// main.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
/* Uncomment the following to detect memery leaks */
//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

#include <iostream>
int main_generate_X_basis(int argc, char** argv);
int main_generate_Reindex(int argc, char** argv);
int main_generate_E4t(int argc, char** argv);
int main_generate_E4_to_E4t(int argc, char** argv);

int main(int argc, char** argv)
{
	int return_code = main_generate_E4t(argc, argv);

	/* Uncomment the following to detect memery leaks */
	//_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_DEBUG);
	//_CrtDumpMemoryLeaks();
	return return_code;
}
