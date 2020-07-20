#include "dababase.h"
#include "mymath.h"
#include <sstream>


struct Generator
{
	int id;
	std::string name;
	int s;
	int t;
	int v;

	Generator(int gen_id, const char* gen_name, int ss, int tt, int vv) : id(gen_id), name(gen_name), s(ss), t(tt), v(vv) {};
};

/* Execute simple commands.*/
void execute_cmd(sqlite3* conn, const char* cmd)
{
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(conn, cmd, strlen(cmd) + 1, &stmt, NULL);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
}

std::vector<int> str_to_mon(const char* str_mon)
{
	std::vector<int> result;
	if (str_mon[0] == '\0')
		return std::vector<int>();
	std::stringstream ss(str_mon);
	while (ss.good()) {
		int i;
		ss >> i;
		result.push_back(i);
		if (ss.peek() == ',')
			ss.ignore();
	}
	return result;
}

std::string mon_to_str(const std::vector<int>& mon)
{
	std::stringstream ss;
	for (size_t i = 0; i < mon.size(); i++) {
		ss << mon[i];
		if (i < mon.size() - 1)
			ss << ",";
	}
	return ss.str();
}

struct BasisMon
{
	std::vector<int> mon;
	int s;
	int t;
	int v;

	BasisMon(const char* str_mon, int ss, int tt, int vv) : mon(str_to_mon(str_mon)), s(ss), t(tt), v(vv) {};
	BasisMon(std::vector<int> mon_, int ss, int tt, int vv) : mon(mon_), s(ss), t(tt), v(vv) {};
};

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

void load_generators(sqlite3* conn, const char* table_name, std::vector<Generator>& generators)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT gen_id, gen_name, s, t, v FROM ") + table_name + ";";
	sqlite3_prepare_v2(conn, cmd.c_str(), cmd.length() + 1, &stmt, NULL);
	while (sqlite3_step(stmt) == SQLITE_ROW)
		generators.emplace_back(sqlite3_column_int(stmt, 0), sqlite3_column_str(stmt, 1), sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4));
	sqlite3_finalize(stmt);
}

void load_leading_terms(sqlite3* conn, const char* table_name, std::vector<std::vector<std::vector<int>>>& leadings)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT leading_term FROM ") + table_name + ";";
	sqlite3_prepare_v2(conn, cmd.c_str(), cmd.length() + 1, &stmt, NULL);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		const char* str_leading_term = sqlite3_column_str(stmt, 0);
		std::vector<int> mon(str_to_mon(str_leading_term));
		if ((unsigned int)mon[0] >= leadings.size())
			leadings.resize(mon[0] + 1);
		leadings[mon[0]].push_back(mon);
	}
	sqlite3_finalize(stmt);
}

void load_basis(sqlite3* conn, const char* table_name, std::vector<BasisMon>& basis)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT mon, s, t, v FROM ") + table_name + ";";
	sqlite3_prepare_v2(conn, cmd.c_str(), cmd.length() + 1, &stmt, NULL);
	while (sqlite3_step(stmt) == SQLITE_ROW)
		basis.emplace_back(sqlite3_column_str(stmt, 0), sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3));
	sqlite3_finalize(stmt);
}

void load_mon_indices(sqlite3* conn, const char* table_name, std::vector<int>& num_mons)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT t, num_mons FROM ") + table_name + " ORDER BY t;";
	sqlite3_prepare_v2(conn, cmd.c_str(), cmd.length() + 1, &stmt, NULL);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		int t = sqlite3_column_int(stmt, 0);
		int index = sqlite3_column_int(stmt, 1);
		int t1 = num_mons.size();
		if (t > t1)
			num_mons.resize(t, t1 > 0 ? num_mons[t1 - 1] : 0);
		num_mons.push_back(index);
	}
	sqlite3_finalize(stmt);
}

void generate_basis(sqlite3* conn, const char* table_name_basis, const std::vector<Generator>& generators, 
	const std::vector<std::vector<std::vector<int>>>& leadings, int t_max)
{
	/* Compile SQL statements */
	sqlite3_stmt* stmt_update_basis;
	std::string cmd_update_basis = std::string("INSERT INTO ") + table_name_basis + " (mon, s, t, v) VALUES (?1, ?2, ?3, ?4);";
	sqlite3_prepare_v2(conn, cmd_update_basis.c_str(), cmd_update_basis.length() + 1, &stmt_update_basis, NULL);

	/* Load the basis and indices. num_mons[t] is the number of monomials in degree <= t. */
	std::vector<int> num_mons;
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
		basis.push_back(BasisMon(std::vector<int>(), 0, 0, 0));
		num_mons.push_back(1);
		t = 1;

		std::string cmd = std::string("INSERT INTO ") + table_name_basis + " (mon_id, mon, s, t, v) VALUES (0, \"\", 0, 0, 0);";
		execute_cmd(conn, cmd.c_str());
	}

	/* Add new basis */
	for (; t <= t_max; t++) {
		std::cout << "t=" << t << '\n';
		for (auto gen : generators) {
			int t1 = t - gen.t;
			if (t1 >= 0) {
				std::cout << "  gen=" << gen.name << '\n';
				size_t index1 = t1 > 0 ? num_mons[t1 - 1] : 0;
				size_t index2 = num_mons[t1];
				for (size_t i = index1; i < index2; i++) {
					if (basis[i].mon.empty() || gen.id <= basis[i].mon[0]) {
						std::vector<int> mon(mul(basis[i].mon, { gen.id, 1 }));
						if ((unsigned int)gen.id >= leadings.size() || std::none_of(leadings[gen.id].begin(), leadings[gen.id].end(),
							[mon](std::vector<int> _m) { return divides(_m, mon); })) {
							int s = basis[i].s + gen.s;
							int v = basis[i].v + gen.v;
							basis.emplace_back(mon, s, t, v);
						}
					}
				}
			}
		}
		std::sort(basis.begin() + num_mons.back(), basis.end(), cmp_basis);
		int num_mons_prev = num_mons.back();
		num_mons.push_back(basis.size());

		/* Insert the new basis in degree t into the database */
		execute_cmd(conn, "BEGIN TRANSACTION");
		for (size_t i = num_mons_prev; i < basis.size(); i++) {
			std::string str_mon(mon_to_str(basis[i].mon));
			sqlite3_bind_text(stmt_update_basis, 1, str_mon.c_str(), str_mon.size() + 1, NULL);
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

	const char *table_generators, *table_relations, *table_basis;
	int t_max;
	if (argc == 1) {
		table_generators = "E2_generators";
		table_relations = "E2_relations";
		table_basis = "E2_basis";
		t_max = 50;
	}
	else {
		table_generators = argv[1];
		table_relations = argv[2];
		table_basis = argv[3];
		t_max = int(argv[4]);
	}
	std::vector<Generator> generators;
	load_generators(conn, table_generators, generators);

	std::vector<std::vector<std::vector<int>>> leadings;
	load_leading_terms(conn, table_relations, leadings);

	generate_basis(conn, table_basis, generators, leadings, 30);

	sqlite3_close(conn);
	return 0;
}