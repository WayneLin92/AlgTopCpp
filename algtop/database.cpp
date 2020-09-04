#include "database.h"


void Database::execute_cmd(const std::string& cmd) const
{
	sqlite3_stmt* stmt;
	sqlite3_prepare_v100(cmd, &stmt);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
}

std::vector<Deg> Database::load_gen_degs(const std::string& table_name) const
{
	std::vector<Deg> gen_degs;
	Statement stmt;
	stmt.init(*this, "SELECT s, t, v FROM " + table_name + " ORDER BY gen_id;");
	while (stmt.step() == SQLITE_ROW)
		gen_degs.emplace_back(stmt.column_int(0), stmt.column_int(1), stmt.column_int(2));
	std::cout << "gen_degs loaded from " << table_name << ", size=" << gen_degs.size() << '\n';
	return gen_degs;
}

array3d Database::load_leading_terms(const std::string& table_name) const
{
	array3d leadings;
	Statement stmt;
	stmt.init(*this, "SELECT leading_term FROM " + table_name + ";");
	while (stmt.step() == SQLITE_ROW) {
		array mon(str_to_array(stmt.column_str(0)));
		if (size_t(mon[0]) >= leadings.size())
			leadings.resize(size_t(mon[0]) + 1);
		leadings[mon[0]].push_back(mon);
	}
	std::cout << "leadings loaded from " << table_name << ", size=" << leadings.size() << '\n';
	return leadings;
}

std::map<Deg, array2d> Database::load_basis(const std::string& table_name) const
{
	std::map<Deg, array2d> basis;
	Statement stmt;
	stmt.init(*this, "SELECT s, t, v, mon FROM " + table_name + " ORDER BY mon_id;");
	while (stmt.step() == SQLITE_ROW) {
		Deg d = { stmt.column_int(0), stmt.column_int(1), stmt.column_int(2) };
		basis[d].push_back(str_to_array(stmt.column_str(3)));
	}
	std::cout << "basis loaded from " << table_name << ", size=" << basis.size() << '\n';
	return basis;
}


void Database::save_basis(const std::string& table_name, const std::map<Deg, array2d>& basis) const
{
	Statement stmt;
	stmt.init(*this, "INSERT INTO " + table_name + " (mon, s, t, v) VALUES (?1, ?2, ?3, ?4);");

	execute_cmd("BEGIN TRANSACTION");
	for (auto& [deg, basis_d] : basis) {
		for (auto& m : basis_d) {
			stmt.bind_str(1, array_to_str(m));
			stmt.bind_int(2, deg.s);
			stmt.bind_int(3, deg.t);
			stmt.bind_int(4, deg.v);
			stmt.step();
			stmt.reset();
		}
	}
	execute_cmd("END TRANSACTION");
}