#include "main.h"

/*********** FUNCTIONS **********/

void execute_cmd(sqlite3* conn, const std::string& cmd)
{
	sqlite3_stmt* stmt;
	sqlite3_prepare_v100(conn, cmd, &stmt);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
}

void load_gen_degs(sqlite3* conn, const std::string& table_name, std::vector<Deg>& gen_degs)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT s, t, v FROM ") + table_name + " ORDER BY gen_id;";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		gen_degs.emplace_back(sqlite3_column_int(stmt, 0), sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2));
	}
	sqlite3_finalize(stmt);
}

void load_leading_terms(sqlite3* conn, const std::string& table_name, array3d& leadings)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT leading_term FROM ") + table_name + ";";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		array mon(str_to_array(sqlite3_column_str(stmt, 0)));
		if (size_t(mon[0]) >= leadings.size())
			leadings.resize(size_t(mon[0]) + 1);
		leadings[mon[0]].push_back(mon);
	}
	sqlite3_finalize(stmt);
}

void load_basis(sqlite3* conn, const std::string& table_name, std::map<Deg, array2d>& basis)
{
	sqlite3_stmt* stmt;
	std::string cmd = "SELECT mon, s, t, v FROM " + table_name + " ORDER BY mon_id;";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		Deg d = { sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3) };
		basis[d].push_back(str_to_array(sqlite3_column_str(stmt, 0)));
	}
	sqlite3_finalize(stmt);
}

void generate_basis(const Database& db, const std::string& table_prefix, int t_max, bool drop_existing=false)
{
	/* load gen_degs, leadings, basis */
	if (drop_existing)
		db.execute_cmd("DELETE FROM " + table_prefix + "_basis;");
	std::vector<Deg> gen_degs = db.load_gen_degs(table_prefix + "_generators");
	array3d leadings = db.load_leading_terms(table_prefix + "_relations");
	std::map<Deg, array2d> basis = db.load_basis(table_prefix + "_basis");

	/* starting t value */
	int t_min;
	if (!basis.empty())
		t_min = basis.rbegin()->first.t + 1;
	else {
		basis[Deg{ 0, 0, 0 }].push_back({}); /* if no monomial present insert the unit */
		t_min = 1;
		db.execute_cmd("INSERT INTO " + table_prefix + "_basis (mon_id, mon, s, t, v) VALUES (0, \"\", 0, 0, 0);");
	}

	/* Add new basis */
	for (int t = t_min; t <= t_max; t++) {
		std::map<Deg, array2d> basis_new;
		std::cout << "t=" << t << "          \r";
		for (int gen_id = (int)gen_degs.size() - 1; gen_id >= 0; --gen_id) {
			int t1 = t - gen_degs[gen_id].t;
			if (t1 >= 0) {
				auto p1 = basis.lower_bound(Deg{ 0, t1, 0 });
				auto p2 = basis.lower_bound(Deg{ 0, t1 + 1, 0 });
				for (auto p = p1; p != p2; ++p) {
					for (const auto& m : p->second) {
						if (m.empty() || gen_id <= m[0]) {
							array mon(mul(m, { gen_id, 1 }));
							if ((size_t)gen_id >= leadings.size() || std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(),
								[&mon](const array& _m) { return divides(_m, mon); }))
								basis_new[p->first + gen_degs[gen_id]].push_back(std::move(mon));
						}
					}
				}
			}
		}
		/* Insert the new basis in degree t into the database */
		db.save_basis(table_prefix + "_basis", basis_new);
		basis.merge(basis_new);
	}
}

int main_generate_basis(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");

	std::string table_prefix;
	int t_max;
	if (argc == 1) {
		table_prefix = "E4b6";
		t_max = 74;
	}
	else {
		table_prefix = argv[1];
		t_max = atoi(argv[2]);
	}

	generate_basis(db, table_prefix, t_max, true);
	return 0;
}