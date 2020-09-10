// main.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
/* Uncomment the following to detect memery leaks */
//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

#include "main.h"
#include "myparser.h"
#include <fstream>
#include <sstream>

int main_generate_basis(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	std::string table_prefix = "E4b6";
	int t_max = 74;
	generate_basis(db, table_prefix, t_max, true);
	return 0;
}

int main_generate_diff(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	std::string table_prefix = "E2";
	generate_mon_diffs(db, table_prefix, 2);
	return 0;
}

int main_generate_ss(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	const char* table_basis = "E2_basis", * table_ss = "E2_ss";
	int t_max = 74;
	generate_ss(db, table_basis, table_ss, 2);
	return 0;
}

int main_generate_next_page(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	std::string table_prefix = "E2", table_H_prefix = "E4";
	generate_next_page(db, table_prefix, table_H_prefix, 2);
	return 0;
}

int main(int argc, char** argv)
{
	int return_code = main_test(argc, argv);

	/* Uncomment the following to detect memery leaks */
	//_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_DEBUG);
	//_CrtDumpMemoryLeaks();
	return return_code;
}