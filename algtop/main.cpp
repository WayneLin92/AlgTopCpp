// tmp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
// Uncomment the following to detect memery leaks
//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

#include "mymath.h"
#include "dababase.h"
#include "myparser.h"
#include <fstream>
#include <sstream>

int main(int argc, char** argv)
{
	int return_code = main_generate_diff(argc, argv);

	// Uncomment the following to detect memory leaks
	//_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_DEBUG);
	//_CrtDumpMemoryLeaks();
	return return_code;
}