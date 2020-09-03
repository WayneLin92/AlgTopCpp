#include "database.h"

struct BasisSS
{
	array2d basis_ind;
	array2d diffs_ind;
	array levels;
};

void load_mon_diffs(sqlite3* conn, const std::string& table_name, std::map<Deg, array2d>& diffs)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT diff, s, t, v FROM ") + table_name + " ;";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		const Deg deg = { sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3) };
		diffs[deg].push_back(str_to_array(sqlite3_column_str(stmt, 0)));
	}
	sqlite3_finalize(stmt);
}

void generate_ss(sqlite3* conn, const std::string& table_name_basis, const std::string& table_ss, int r)
{
	/* Compile SQL statements */
	sqlite3_stmt* stmt_update_basis;
	std::string cmd_update_basis = std::string("INSERT INTO ") + table_ss + " (base, diff, level, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5, ?6);";
	sqlite3_prepare_v100(conn, cmd_update_basis, &stmt_update_basis);

	execute_cmd(conn, "DELETE FROM " + table_ss + ";");

	std::map<Deg, array2d> mon_diffs;
	load_mon_diffs(conn, table_name_basis, mon_diffs);
	std::cout << "mon_diffs loaded! Size=" << mon_diffs.size() << '\n';

	std::map<Deg, BasisSS> basis_ss;

	/* fill basis_ss */
	int prev_t = 0;
	for (auto& [deg, mon_diffs_d] : mon_diffs) {
		if (deg.t > prev_t) {
			prev_t = deg.t;
			std::cout << "Compute basis_ss, t=" << deg.t << "          \r";
		}
		BasisSS& basis_ss_d = basis_ss[deg];

		array range_ = range((int)mon_diffs_d.size());

		array lead_image;
		for (const auto& base : basis_ss_d.basis_ind)
			lead_image.push_back(base.front());
		std::sort(lead_image.begin(), lead_image.end());
		if (!lead_image.empty() && (range_.empty() || lead_image.back() > range_.back()))
			std::cout << "test\n";

		array x = add_vectors(range_, lead_image);

		array2d fx;
		for (int xi : x)
			fx.push_back(mon_diffs_d[xi]);

		array2d image, kernel, g;
		set_linear_map(x, fx, image, kernel, g);

		/* fill with other cycles after the boundaries */
		array lead_kernel;
		for (auto& cycle : kernel) {
			lead_kernel.push_back(cycle.front());
			basis_ss_d.basis_ind.push_back(cycle);
			basis_ss_d.diffs_ind.push_back({});
			basis_ss_d.levels.push_back(T_MAX - r);
		}
		std::sort(lead_kernel.begin(), lead_kernel.end());

		/* fill with boundaries */
		Deg deg_diff = deg + Deg{ 1, 0, -r };
		for (auto& boundary : image) {
			basis_ss[deg_diff].basis_ind.push_back(std::move(boundary));
			basis_ss[deg_diff].diffs_ind.push_back({});
			basis_ss[deg_diff].levels.push_back(r);
		}

		/* fill with the rest */
		array rest = add_vectors(x, lead_kernel);
		for (int i : rest) {
			basis_ss_d.basis_ind.push_back({ i });
			basis_ss_d.diffs_ind.push_back(std::move(mon_diffs[deg][i]));
			basis_ss_d.levels.push_back(T_MAX);
		}
	}

	/* Insert into the database */
	execute_cmd(conn, "BEGIN TRANSACTION");
	for (const auto& [deg, basis_ss_d] : basis_ss) {
		for (size_t i = 0; i < basis_ss_d.basis_ind.size(); ++i) {
			sqlite3_bind_str(stmt_update_basis, 1, array_to_str(basis_ss_d.basis_ind[i]));
			sqlite3_bind_str(stmt_update_basis, 2, array_to_str(basis_ss_d.diffs_ind[i]));
			sqlite3_bind_int(stmt_update_basis, 3, basis_ss_d.levels[i]);
			sqlite3_bind_int(stmt_update_basis, 4, deg.s);
			sqlite3_bind_int(stmt_update_basis, 5, deg.t);
			sqlite3_bind_int(stmt_update_basis, 6, deg.v);
			sqlite3_step(stmt_update_basis);
			sqlite3_reset(stmt_update_basis);
		}
	}
	execute_cmd(conn, "END TRANSACTION");

	sqlite3_finalize(stmt_update_basis);
}

int main_generate_ss(int argc, char** argv)
{
	sqlite3* conn;
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)", &conn);

	const char* table_basis, * table_ss;
	int t_max;
	if (argc == 1) {
		table_basis = "E2_basis";
		table_ss = "E2_ss";
		t_max = 74;
	}
	else {
		table_basis = argv[1];
		table_ss = argv[2];
		t_max = atoi(argv[3]);
	}

	generate_ss(conn, table_basis, table_ss, 2);

	sqlite3_close(conn);
	return 0;
}