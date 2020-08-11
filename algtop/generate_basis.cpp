#include "database.h"

/*********** STRUCTS **********/
struct Generator
{
	int id;
	int s;
	int t;
	int v;

	Generator(int gen_id, int s_, int t_, int v_) : id(gen_id), s(s_), t(t_), v(v_) {};
};

struct BasisMon
{
	array mon;
	int s;
	int t;
	int v;

	BasisMon(const char* str_mon, int ss, int tt, int vv) : mon(str_to_array(str_mon)), s(ss), t(tt), v(vv) {};
	BasisMon(array&& mon_, int ss, int tt, int vv) : mon(mon_), s(ss), t(tt), v(vv) {};
};

/*********** FUNCTIONS **********/
/* Execute simple commands.*/
void execute_cmd(sqlite3* conn, const char* cmd)
{
	sqlite3_stmt* stmt;
	sqlite3_prepare_v100(conn, cmd, &stmt);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
}

/* Sort monomials by their s, v and mon by reversed lexicographical order */
bool cmp_basis(const BasisMon& m1, const BasisMon& m2)
{
	if (m1.t < m2.t)
		return true;
	else if (m1.t == m2.t)
		if (m1.s < m2.s)
			return true;
		else if (m1.s == m2.s)
			if (m1.v < m2.v)
				return true;
			else if (m1.v == m2.v)
				return cmp_mons(m1.mon, m2.mon);
	return false;
}

void load_generators(sqlite3* conn, const std::string& table_name, std::vector<Generator>& generators)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT gen_id, s, t, v FROM ") + table_name + ";";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW)
		generators.emplace_back(sqlite3_column_int(stmt, 0), sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3));
	sqlite3_finalize(stmt);
}

void load_leading_terms(sqlite3* conn, const std::string& table_name, std::vector<std::vector<array>>& leadings)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT leading_term FROM ") + table_name + ";";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		const char* str_leading_term = sqlite3_column_str(stmt, 0);
		array mon(str_to_array(str_leading_term));
		if (size_t(mon[0]) >= leadings.size())
			leadings.resize(size_t(mon[0]) + 1);
		leadings[mon[0]].push_back(mon);
	}
	sqlite3_finalize(stmt);
}

void load_basis(sqlite3* conn, const std::string& table_name, std::vector<BasisMon>& basis)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT mon, s, t, v FROM ") + table_name + ";";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW)
		basis.emplace_back(sqlite3_column_str(stmt, 0), sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3));
	sqlite3_finalize(stmt);
}

void load_mon_indices(sqlite3* conn, const char* table_name, array& num_mons)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT t, num_mons FROM ") + table_name + " ORDER BY t;";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		int t = sqlite3_column_int(stmt, 0);
		int index = sqlite3_column_int(stmt, 1);
		size_t t1 = num_mons.size();
		if (t > t1)
			num_mons.resize(t, t1 > 0 ? num_mons[t1 - 1] : 0);
		num_mons.push_back(index);
	}
	sqlite3_finalize(stmt);
}

void generate_basis(sqlite3* conn, const std::string& table_name_basis, const std::vector<Generator>& generators, 
	const array3d& leadings, int t_max)
{
	/* Compile SQL statements */
	sqlite3_stmt* stmt_update_basis;
	std::string cmd_update_basis = std::string("INSERT INTO ") + table_name_basis + " (mon, s, t, v) VALUES (?1, ?2, ?3, ?4);";
	sqlite3_prepare_v100(conn, cmd_update_basis, &stmt_update_basis);

	/* Load the basis and indices. num_mons[t] is the number of monomials in degree <= t. */
	array num_mons;
	std::string table_indices_t = std::string(table_name_basis) + "_indices_t";
	load_mon_indices(conn, table_indices_t.c_str(), num_mons);
	std::vector<BasisMon> basis;
	load_basis(conn, table_name_basis, basis);

	/* Get the current largest t value */
	int t;
	if (!basis.empty())
		t = basis[basis.size() - 1].t + 1;
	else {
		/* If no monomial present insert the unit.*/
		basis.push_back(BasisMon(array(), 0, 0, 0));
		num_mons.push_back(1);
		t = 1;
		std::string cmd = std::string("INSERT INTO ") + table_name_basis + " (mon_id, mon, s, t, v) VALUES (0, \"\", 0, 0, 0);";
		execute_cmd(conn, cmd.c_str());
	}

	/* Add new basis */
	for (; t <= t_max; t++) {
		std::cout << "t=" << t << '\n';
		for (const auto& gen : generators) {
			int t1 = t - gen.t;
			if (t1 >= 0) {
				size_t index1 = t1 > 0 ? num_mons[size_t(t1) - 1] : 0;
				size_t index2 = num_mons[t1];
				for (size_t i = index1; i < index2; i++) {
					if (basis[i].mon.empty() || gen.id <= basis[i].mon[0]) {
						array mon(mul(basis[i].mon, { gen.id, 1 }));
						if ((size_t)gen.id >= leadings.size() || std::none_of(leadings[gen.id].begin(), leadings[gen.id].end(),
							[mon](array _m) { return divides(_m, mon); })) {
							int s = basis[i].s + gen.s;
							int v = basis[i].v + gen.v;
							basis.emplace_back(std::move(mon), s, t, v);
						}
					}
				}
			}
		}
		int num_mons_prev = num_mons.back();
		std::sort(basis.begin() + num_mons_prev, basis.end(), cmp_basis);
		num_mons.push_back(int(basis.size()));

		/* Insert the new basis in degree t into the database */
		execute_cmd(conn, "BEGIN TRANSACTION");
		for (size_t i = num_mons_prev; i < basis.size(); i++) {
			sqlite3_bind_str(stmt_update_basis, 1, array_to_str(basis[i].mon));
			sqlite3_bind_int(stmt_update_basis, 2, basis[i].s);
			sqlite3_bind_int(stmt_update_basis, 3, basis[i].t);
			sqlite3_bind_int(stmt_update_basis, 4, basis[i].v);
			sqlite3_step(stmt_update_basis);
			sqlite3_reset(stmt_update_basis);
		}
		execute_cmd(conn, "END TRANSACTION");
	}
	sqlite3_finalize(stmt_update_basis);
}

int main_generate_basis(int argc, char** argv)
{
	sqlite3* conn;
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\ss.db)", &conn);

	std::string table_prefix;
	int t_max;
	if (argc == 1) {
		table_prefix = "E4";
		t_max = 200;
	}
	else {
		table_prefix = argv[1];
		t_max = atoi(argv[2]);
	}
	std::vector<Generator> generators;
	load_generators(conn, table_prefix + "_generators", generators);

	std::vector<std::vector<array>> leadings;
	load_leading_terms(conn, table_prefix + "_relations", leadings);

	generate_basis(conn, table_prefix + "_basis", generators, leadings, t_max);

	sqlite3_close(conn);
	return 0;
}