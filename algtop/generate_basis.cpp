#include "database.h"

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
	std::string cmd = std::string("SELECT mon, s, t, v FROM ") + table_name + " ORDER BY mon_id;";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		Deg d = { sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3) };
		basis[d].push_back(str_to_array(sqlite3_column_str(stmt, 0)));
	}
	sqlite3_finalize(stmt);
}

void generate_basis(sqlite3* conn, const std::string& table_prefix, int t_max, bool drop_existing=false)
{
	/* load gen_degs, leadings, basis */

	std::vector<Deg> gen_degs;
	load_gen_degs(conn, table_prefix + "_generators", gen_degs);
	std::cout << "gen_degs loaded! Size=" << gen_degs.size() << '\n';

	array3d leadings;
	load_leading_terms(conn, table_prefix + "_relations", leadings);
	std::cout << "leadings loaded! Size=" << leadings.size() << '\n';

	std::map<Deg, array2d> basis;
	if (drop_existing)
		execute_cmd(conn, "DELETE FROM " + table_prefix + "_basis;");
	else {
		load_basis(conn, table_prefix + "_basis", basis);
		std::cout << "basis loaded! Size=" << basis.size() << '\n';
	}

	/* Get the starting t value */

	int t_min;
	if (!basis.empty())
		t_min = basis.rbegin()->first.t + 1;
	else {
		basis[Deg{ 0, 0, 0 }].push_back({}); /* if no monomial present insert the unit */
		t_min = 1;
		std::string cmd = "INSERT INTO " + table_prefix + "_basis (mon_id, mon, s, t, v) VALUES (0, \"\", 0, 0, 0);";
		execute_cmd(conn, cmd.c_str());
	}

	/* Add new basis */

	sqlite3_stmt* stmt_update_basis;
	std::string cmd_update_basis = "INSERT INTO " + table_prefix + "_basis (mon, s, t, v) VALUES (?1, ?2, ?3, ?4);";
	sqlite3_prepare_v100(conn, cmd_update_basis, &stmt_update_basis);

	for (int t = t_min; t <= t_max; t++) {
		std::map<Deg, array2d> basis_new;
		std::cout << "t=" << t << "          \r";
		for (int i = 0; i < (int)gen_degs.size(); ++i) {
			int t1 = t - gen_degs[i].t;
			if (t1 >= 0) {
				auto p1 = basis.lower_bound(Deg{ 0, t1, 0 });
				auto p2 = basis.lower_bound(Deg{ 0, t1 + 1, 0 });
				for (auto p = p1; p != p2; ++p) {
					for (const auto& m : p->second) {
						if (m.empty() || i <= m[0]) {
							array mon(mul(m, { i, 1 }));
							if ((size_t)i >= leadings.size() || std::none_of(leadings[i].begin(), leadings[i].end(),
								[&mon](array _m) { return divides(_m, mon); })) {
								int s = p->first.s + gen_degs[i].s;
								int v = p->first.v + gen_degs[i].v;
								basis_new[Deg{ s, t, v }].push_back(std::move(mon));
							}
						}
					}
				}
			}
		}

		/* Insert the new basis in degree t into the database */

		execute_cmd(conn, "BEGIN TRANSACTION");
		for (auto p = basis_new.begin(); p != basis_new.end(); ++p) {
			for (const auto& m : p->second) {
				sqlite3_bind_str(stmt_update_basis, 1, array_to_str(m));
				sqlite3_bind_int(stmt_update_basis, 2, p->first.s);
				sqlite3_bind_int(stmt_update_basis, 3, p->first.t);
				sqlite3_bind_int(stmt_update_basis, 4, p->first.v);
				sqlite3_step(stmt_update_basis);
				sqlite3_reset(stmt_update_basis);
			}
		}
		execute_cmd(conn, "END TRANSACTION");
		basis.merge(basis_new);
	}
	sqlite3_finalize(stmt_update_basis);
}

int main_generate_basis(int argc, char** argv)
{
	sqlite3* conn;
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)", &conn);

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

	generate_basis(conn, table_prefix, t_max, true);

	sqlite3_close(conn);
	return 0;
}