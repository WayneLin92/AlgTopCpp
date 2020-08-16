#include "database.h"

/*
** Load E4_relations
** compute Ann(h_0)
*/

void load_gb(sqlite3* conn, const std::string& table_name, array3d& gb)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT leading_term, basis FROM ") + table_name + ";";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		const char* str_leading_term = sqlite3_column_str(stmt, 0);
		const char* str_basis = sqlite3_column_str(stmt, 1);
		array lead(str_to_array(str_leading_term));
		array2d basis(str_to_array2d(str_basis));
		array2d g = add(basis, { lead });
		gb.push_back(std::move(g));
	}
	sqlite3_finalize(stmt);
}

void load_degs(sqlite3* conn, const std::string& table_name, array& gen_degs)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT t FROM ") + table_name + ";";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		gen_degs.push_back(sqlite3_column_int(stmt, 0));
	}
	sqlite3_finalize(stmt);
}

void generate_E4t(sqlite3* conn, const std::string& table_prefix, const std::string& table_H_prefix, int r)
{

}

int main_generate_E4t(int argc, char** argv)
{
	sqlite3* conn;
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\ss.db)", &conn);

	array gen_degs;
	load_degs(conn, "E4_generators", gen_degs);
	std::cout << gen_degs.size() << '\n';

	array3d gb;
	load_gb(conn, "E4_relations", gb);
	std::cout << gb.size() << '\n';

	array3d polys = { {{0, 1}} };
	array4d result = ann_seq(gb, polys, gen_degs, 200);
	std::cout << result << '\n';
	std::cout << result.size() << '\n';

	array dims;
	dims.resize(201);
	for (const array3d& v : result) {
		dims[deg(v[0], gen_degs)]++;
		if (deg(v[0], gen_degs) == 20)
			std::cout << v[0] << '\n';
	}
	/*for (int i = 0; i < dims.size(); i++) {
		std::cout << i << ": " << dims[i] << '\n';
	}*/

	return 0;
}