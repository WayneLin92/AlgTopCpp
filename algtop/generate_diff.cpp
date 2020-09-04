#include "main.h"

/********** STRUCTS AND CLASSES **********/

struct DgaBasis
{
	array2d basis;
	array diffs_determined;
	array2d diffs_indices;
	array3d diffs;
};

/********** FUNCTIONS **********/

void load_gen_diffs(sqlite3* conn, const std::string& table_name, array3d& diffs)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT gen_diff FROM ") + table_name + ";";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW)
		diffs.push_back(str_to_array2d(sqlite3_column_str(stmt, 0)));
	sqlite3_finalize(stmt);
}

void load_gb(sqlite3* conn, const std::string& table_name, array3d& gb)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT leading_term, basis FROM ") + table_name + ";";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		array lead(str_to_array(sqlite3_column_str(stmt, 0)));
		array2d basis(str_to_array2d(sqlite3_column_str(stmt, 1)));
		array2d g = add(basis, { lead });
		gb.push_back(std::move(g));
	}
	sqlite3_finalize(stmt);
}

void load_dga_basis(sqlite3* conn, const std::string& table_name, std::map<Deg, DgaBasis>& basis, int r)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT mon, diff, s, t, v FROM ") + table_name + " ORDER BY mon_id;";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	int prev_t = 0;
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		Deg d = { sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4) };
		if (d.t > prev_t) {
			std::cout << "loading_dga_basis, t=" << d.t << "          \r";
			prev_t = d.t;
		}
		if (sqlite3_column_type(stmt, 1) == SQLITE_NULL) {
			basis[d].basis.push_back(str_to_array(sqlite3_column_str(stmt, 0)));
			basis[d].diffs_determined.push_back(false);
			basis[d].diffs_indices.push_back(array());
			basis[d].diffs.push_back(array2d());
		}
		else {
			basis[d].basis.push_back(str_to_array(sqlite3_column_str(stmt, 0)));
			basis[d].diffs_determined.push_back(true);
			basis[d].diffs_indices.push_back(str_to_array(sqlite3_column_str(stmt, 1)));
			basis[d].diffs.push_back(array2d());
		}
	}
	for (auto& [d, basis_d] : basis) {
		for (size_t i = 0; i < basis_d.basis.size(); ++i)
			if (basis_d.diffs_determined[i])
				basis_d.diffs[i] = indices_to_poly(basis_d.diffs_indices[i], basis[d + Deg{ 1, 0, -r }].basis);
	}
	sqlite3_finalize(stmt);
}

void generate_mon_diffs(sqlite3* conn, const std::string& table_prefix, int r)
{
	/* load gen_degs, gen_diffs, gb and basis */

	std::vector<Deg> gen_degs;
	load_gen_degs(conn, table_prefix + "_generators", gen_degs);
	std::cout << "gen_degs loaded! Size=" << gen_degs.size() << '\n';

	array3d gen_diffs;
	load_gen_diffs(conn, table_prefix + "_generators", gen_diffs);
	std::cout << "gen_diffs loaded! Size=" << gen_diffs.size() << '\n';

	array3d gb;
	load_gb(conn, table_prefix + "_relations", gb);
	std::cout << "gb loaded! Size=" << gb.size() << '\n';

	std::map<Deg, DgaBasis> basis;
	load_dga_basis(conn, table_prefix + "_basis", basis, r);
	std::cout << "basis loaded! Size=" << basis.size() << '\n';

	/* compute diffs */

	sqlite3_stmt* stmt_update_basis;
	std::string cmd_update_basis = std::string("UPDATE ") + table_prefix + "_basis SET diff=?1 WHERE mon=?2;";
	sqlite3_prepare_v100(conn, cmd_update_basis, &stmt_update_basis);

	execute_cmd(conn, "BEGIN TRANSACTION");
	int prev_t = 0;
	for (auto& [d, basis_d] : basis) {
		for (size_t i = 0; i < basis_d.basis.size(); ++i) {
			if (!basis_d.diffs_determined[i]) {
				array diff_indices;
				if (d.t > prev_t) {
					prev_t = d.t;
					std::cout << "computing diffs, t=" << d.t << "          \r";
					execute_cmd(conn, "END TRANSACTION");
					execute_cmd(conn, "BEGIN TRANSACTION");
				}
				if (basis_d.basis[i].empty()) {
					basis_d.diffs[i] = array2d();
					basis_d.diffs_determined[i] = true;
				}
				else {
					int gen_id = basis_d.basis[i].front();
					array mon1 = div(basis_d.basis[i], { gen_id, 1 });
					Deg d1 = d - gen_degs[gen_id];
					size_t index_mon1 = std::lower_bound(basis.at(d1).basis.begin(), basis.at(d1).basis.end(), mon1, cmp_mons) - basis.at(d1).basis.begin();
					basis_d.diffs[i] = reduce(add(mul(gen_diffs[gen_id], mon1), mul({ { gen_id, 1 } }, basis.at(d1).diffs[index_mon1])), gb);
					Deg d_diff = d + Deg{ 1, 0, -r };
					diff_indices = poly_to_indices(basis_d.diffs[i], basis[d_diff].basis);
					basis_d.diffs_determined[i] = true;
				}
				std::string str_diff(array_to_str(diff_indices));
				sqlite3_bind_str(stmt_update_basis, 1, str_diff);
				sqlite3_bind_str(stmt_update_basis, 2, array_to_str(basis_d.basis[i]));
				sqlite3_step(stmt_update_basis);
				sqlite3_reset(stmt_update_basis);
			}
		}
	}
	execute_cmd(conn, "END TRANSACTION");
	sqlite3_finalize(stmt_update_basis);
}

int test1()
{
	return 0;
}

int main_generate_diff(int argc, char** argv)
{
	//return test1();

	sqlite3* conn;
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)", &conn);

	std::string table_prefix;
	if (argc == 1) {
		table_prefix = "E2";
	}
	else {
		table_prefix = argv[1];
	}

	generate_mon_diffs(conn, table_prefix, 2);

	sqlite3_close(conn);
	return 0;
}