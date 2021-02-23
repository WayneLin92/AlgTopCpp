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
	if (argc == 6) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		const char* table_prefix = argv[3];
		int t_max = std::stoi(argv[4]);
 		bool drop_existing = std::stoi(argv[5]) == 0 ? false : true;

		Database db(db_path.c_str());
		generate_basis(db, table_prefix, t_max, drop_existing);
		return 0;
	}
	else {
		std::cout << "usage ss basis <db_path> <table_name> <t_max> <drop_existing : 1 | 0>\n";
		return 0;
	}
}

int main_generate_diff(int argc, char** argv)
{
	if (argc == 5) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		const char* table_prefix = argv[3];
		int r = std::stoi(argv[4]);

		Database db(db_path.c_str());
		generate_mon_diffs(db, table_prefix, r);
		return 0;
	}
	else {
		std::cout << "usage ss diff <db_path> <table_name> <r>\n";
		return 0;
	}
}

int main_generate_ss(int argc, char** argv)
{
	if (argc == 5) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		std::string table_prefix = argv[3];
		int r = std::stoi(argv[4]);

		Database db(db_path.c_str());
		generate_ss(db, table_prefix + "_basis", table_prefix + "_ss", r);
		return 0;
	}
	else {
		std::cout << "usage ss ss <db_path> <table_name> <r>\n";
		return 0;
	}
}

int main_generate_next_page(int argc, char** argv)
{
	if (argc == 6) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		std::string table_prefix = argv[3];
		std::string table_H_prefix = argv[4];
		int r = std::stoi(argv[5]);

		Database db(db_path.c_str());
		generate_next_page(db, table_prefix, table_H_prefix, r);
		return 0;
	}
	else {
		std::cout << "usage: ss next_page <db_path> <table_name> <tables_name_H> <r>\n";
		return 0;
	}
}

int main_deduce_zero_diffs(int argc, char** argv)
{
	if (argc == 5) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		std::string table_prefix = argv[3];
		int t_max = std::stoi(argv[4]);

		Database db(db_path.c_str());
		DiffType d_type = GetDiffType(table_prefix);
		WrapDeduceZeroDiffs(db, table_prefix, t_max);
		return 0;
	}
	else {
		std::cout << "usage: ss deduce_zero_diffs <db_path> <table_name> <t_max> [<bForEt> =0]\n";
		return 0;
	}
}

int main_deduce_diffs_for_Et(int argc, char** argv)
{
	if (argc == 5) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		std::string table_prefix = argv[3];
		int t_max = std::stoi(argv[4]);

		Database db(db_path.c_str());
		WrapDeduceDiffsForEt(db, table_prefix, t_max);
		return 0;
	}
	else {
		std::cout << "usage: ss deduce_for_Et <db_path> <table_name> <t_max>\n";
		return 0;
	}
}

int main_deduce_diffs_by_trying(int argc, char** argv)
{
	if (argc == 6) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		std::string table_prefix = argv[3];
		int t_try = std::stoi(argv[4]);
		int t_test = std::stoi(argv[5]);

		Database db(db_path.c_str());
		WrapDeduceDiffsByTrying(db, table_prefix, t_try, t_test);
		return 0;
	}
	else {
		std::cout << "usage: ss deduce_by_trying <db_path> <table_name> <t_try> <t_test> [<iter_max> =100 <bForEt> =0]\n";
		return 0;
	}
}

int main_deduce_diffs_by_naturality(int argc, char** argv)
{
	if (argc == 8) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		std::string table_prefix = argv[3];
		std::string table_prefix1 = argv[4];
		std::string column = argv[5];
		int t_try = std::stoi(argv[6]);
		int t_test = std::stoi(argv[7]);

		Database db(db_path.c_str());
		WrapDeduceDiffsByNat(db, table_prefix, table_prefix1, column, t_try, t_test);
		return 0;
	}
	else {
		std::cout << "usage: ss deduce_by_nat <db_path> <table_name> <table_name1> <column_image> <t_try> <t_test> [<iter_max> =100 <bForEt> =0]\n";
		return 0;
	}
}

int main_check_diffs_by_naturality(int argc, char** argv)
{
	if (argc == 8) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		std::string table_prefix = argv[3];
		std::string table_prefix1 = argv[4];
		std::string column = argv[5];
		int t_try = std::stoi(argv[6]);
		int t_test = std::stoi(argv[7]);

		Database db(db_path.c_str());
		WrapCheckDiffsByNat(db, table_prefix, table_prefix1, column, t_try, t_test);
		return 0;
	}
	else {
		std::cout << "usage: ss deduce_by_nat <db_path> <table_name> <table_name1> <column_image> <t_try> <t_test> [<iter_max> =100 <bForEt> =0]\n";
		return 0;
	}
}

int main_reset_E4_ss(int argc, char** argv)
{
	Database db("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\tmp.db");
	db.execute_cmd("DELETE FROM E4_ss;");
	db.execute_cmd("INSERT INTO E4_ss SELECT * FROM E4_ss_backup");
	return 0;
}

int main_reset_E4t1_ss(int argc, char** argv)
{
	Database db("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\tmp.db");
	db.execute_cmd("DELETE FROM E4t1_ss;");
	db.execute_cmd("INSERT INTO E4t1_ss SELECT * FROM E4t_ss_tmp");
	return 0;
}

int test_deduce_zero_diffs()
{
	char argv0[100] = "";
	char argv1[100] = "";
	char argv2[100] = "tmp.db";
	char table[100] = "E4";
	char t_max[100] = "230";
	char* argv[] = { argv0, argv1, argv2, table, t_max};
	main_deduce_zero_diffs(5, argv);
	return 0;
}

int test_deduce_diffs_for_Et()
{
	char argv0[100] = "";
	char argv1[100] = "";
	char argv2[100] = "tmp.db";
	char table[100] = "E4h0";
	char t_max[100] = "460";
	char* argv[] = { argv0, argv1, argv2, table, t_max };
	main_deduce_diffs_for_Et(5, argv);
 	return 0;
}

int test_deduce_diffs_by_trying()
{
	char argv0[100] = "";
	char argv1[100] = "";
	char argv2[100] = "tmp.db";
	char table[100] = "E4h0";
	char t_try[100] = "460";
	char t_test[100] = "460";
	char* argv[] = { argv0, argv1, argv2, table, t_try, t_test };
	main_deduce_diffs_by_trying(6, argv);
	return 0;
}

int test_deduce_diffs_by_nat()
{
	char argv0[100] = "";
	char argv1[100] = "";
	char argv2[100] = "tmp.db";
	char argv3[100] = "E4";
	char argv4[100] = "E4h0";
	char argv5[100] = "loc";
	char t_try[100] = "230";
	char t_test[100] = "230";
	char* argv[] = { argv0, argv1, argv2, argv3, argv4, argv5, t_try, t_test};
	main_deduce_diffs_by_naturality(8, argv);
	return 0;
}

int test_check_diffs_by_nat()
{
	char argv0[100] = "";
	char argv1[100] = "";
	char argv2[100] = "tmp.db";
	char argv3[100] = "E4";
	char argv4[100] = "E4t";
	char argv5[100] = "to_E4t";
	char t_try[100] = "189";
	char t_test[100] = "189";
	char* argv[] = { argv0, argv1, argv2, argv3, argv4, argv5, t_try, t_test };
	main_check_diffs_by_naturality(8, argv);
	return 0;
}

int test_deduce_diffs()
{
	Database db("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\tmp.db");
	std::string table_prefix = "E4t1";
	int t_max = 50;
	DeduceDiffs(db, table_prefix, t_max);
	return 0;
}

int test_generate_ss()
{
	char argv0[100] = "";
	char argv1[100] = "";
	char argv2[100] = "tmp.db";
	char argv3[100] = "E4t";
	char argv4[100] = "4";
	char* argv[] = { argv0, argv1, argv2, argv3, argv4 };
	main_generate_ss(5, argv);

	return 0;
}

int test_generate_next_page()
{
	char argv0[100] = "";
	char argv1[100] = "";
	char argv2[100] = "tmp.db";
	char argv3[100] = "E2";
	char argv4[100] = "E4";
	char argv5[100] = "2";
	char* argv[] = { argv0, argv1, argv2, argv3, argv4, argv5 };
	main_generate_next_page(6, argv);

	return 0;
}

int main(int argc, char** argv)
{
	//return test_deduce_zero_diffs(); ////
	return test_deduce_diffs_by_nat(); ////

	if (argc > 1) {
		if (strcmp(argv[1], "basis") == 0)
			return main_generate_basis(argc, argv);
		else if (strcmp(argv[1], "diff") == 0)
			return main_generate_diff(argc, argv);
		else if (strcmp(argv[1], "ss") == 0)
			return main_generate_ss(argc, argv);
		else if (strcmp(argv[1], "next_page") == 0)
			return main_generate_next_page(argc, argv);
		else if (strcmp(argv[1], "deduce_zero_diffs") == 0)
			return main_deduce_zero_diffs(argc, argv);
		else if (strcmp(argv[1], "deduce_for_Et") == 0)
			return main_deduce_diffs_for_Et(argc, argv);
		else if (strcmp(argv[1], "deduce_by_trying") == 0)
			return main_deduce_diffs_by_trying(argc, argv);
		else if (strcmp(argv[1], "reset_E4_ss") == 0)
			return main_reset_E4_ss(argc, argv);
	}
	std::cout << "usage ss <basis, diff, ss, next_page, deduce_zero_diffs, deduce_for_Et, deduce_by_trying, reset_E4_ss> [<args>]\n";
	return 0;

	/* Uncomment the following to detect memery leaks */
	//_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_DEBUG);
	//_CrtDumpMemoryLeaks();
}
