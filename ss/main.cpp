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

int main_deduce_diffs(int argc, char** argv)
{
	if (argc == 9) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		std::string table_prefix = argv[3];
		int t_try = std::stoi(argv[4]);
		int t_test = std::stoi(argv[5]);
		int max_poss = std::stoi(argv[6]);
		int nCheck = std::stoi(argv[7]);
		int seconds = std::stoi(argv[8]);

		Database db(db_path.c_str());
		WrapDeduceDiffs(db, table_prefix, t_try, t_test, max_poss, nCheck, (double)seconds);
		return 0;
	}
	else {
		std::cout << "usage: ss deduce_diffs <db_path> <table_name> <t_try> <t_test> <max_poss> <nCheck> <seconds>\n";
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
		std::cout << "usage: ss deduce_by_nat <db_path> <table_name> <table_name1> <column_image> <t_try> <t_test>\n";
		return 0;
	}
}

int main_check_deduce_diffs_by_naturality_v2(int argc, char** argv)
{
	if (argc == 7) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		std::string table_prefix = argv[3];
		std::string table_prefix1 = argv[4];
		std::string column = argv[5];
		int t_max = std::stoi(argv[6]);

		Database db(db_path.c_str());
		WrapDeduceDiffsByNatV2(db, table_prefix, table_prefix1, column, t_max);
		return 0;
	}
	else {
		std::cout << "usage: ss deduce_by_nat_v2 <db_path> <table_name> <table_name1> <column_image> <t_max>\n";
		return 0;
	}
}

int main_deduce_diffs_for_E4(int argc, char** argv)
{
	if (argc == 5) {
		auto db_path = std::string("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\") + argv[2];
		int t_try = std::stoi(argv[3]);
		int t_test = std::stoi(argv[4]);

		Database db(db_path.c_str());
		WrapDeduceDiffsForE4(db, t_try, t_test);
		return 0;
	}
	else {
		std::cout << "usage: ss deduce_for_E4 <db_path> <t_try> <t_test>\n";
		return 0;
	}
}

int main_reset_tmp1_E4t_ss(int argc, char** argv)
{
	Database db("C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\database\\tmp1.db");
	db.execute_cmd("DELETE FROM E4t_ss;");
	db.execute_cmd("INSERT INTO E4t_ss SELECT * FROM E4t_ss_backup");
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
	char table[100] = "E4t";
	char t_try[100] = "60";
	char t_test[100] = "100";
	char nCheck[100] = "100";
	char* argv[] = { argv0, argv1, argv2, table, t_try, t_test, nCheck };
	main_deduce_diffs(7, argv);
	return 0;
}

int test_deduce_diffs_by_trying1()
{
	char argv0[100] = "";
	char argv1[100] = "";
	char argv2[100] = "tmp.db";
	char table[100] = "E4";
	char t_try[100] = "150";
	char t_test[100] = "230";
	char nCheck[100] = "100";
	char* argv[] = { argv0, argv1, argv2, table, t_try, t_test, nCheck };
	main_deduce_diffs(7, argv);
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

int test_deduce_diffs_by_nat_v2()
{
	char argv0[100] = "";
	char argv1[100] = "";
	char argv2[100] = "tmp.db";
	char argv3[100] = "E4";
	char argv4[100] = "E4t";
	char argv5[100] = "to_E4t";
	char t_max[100] = "189";
	char* argv[] = { argv0, argv1, argv2, argv3, argv4, argv5, t_max };
	main_check_deduce_diffs_by_naturality_v2(7, argv);
	return 0;
}

int test_deduce_diffs_for_E4()
{
	char argv0[100] = "";
	char argv1[100] = "";
	char argv2[100] = "tmp.db";
	char t_try[100] = "230";
	char t_test[100] = "230";
	char* argv[] = { argv0, argv1, argv2, t_try, t_test };
	main_deduce_diffs_for_E4(5, argv);
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
	return test_deduce_diffs_by_trying(); ////
	//return test_deduce_diffs_for_E4(); ////

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
		else if (strcmp(argv[1], "deduce_diffs") == 0)
			return main_deduce_diffs(argc, argv);
		else if (strcmp(argv[1], "reset_E4_ss") == 0)
			return main_reset_tmp1_E4t_ss(argc, argv);
	}
	std::cout << "usage ss <basis, diff, ss, next_page, deduce_zero_diffs, deduce_for_Et, deduce_diffs, reset_E4_ss> [<args>]\n";
	return 0;

	/* Uncomment the following to detect memery leaks */
	//_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_DEBUG);
	//_CrtDumpMemoryLeaks();
}
