#include "dababase.h"
#include "mymath.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <array>


struct DgaGenerator
{
	int id;
	std::string name;
	PolyRaw diff;
	int s;
	int t;
	int v;

	DgaGenerator(int gen_id, const char* gen_name, PolyRaw&& gen_diff, int ss, int tt, int vv) :
		id(gen_id), name(gen_name), diff(gen_diff), s(ss), t(tt), v(vv) {};
};

struct DgaBasisMon
{
	MonRaw mon;
	bool diff_determined;
	PolyRaw diff;
	std::vector<int> diff_indices;
	int s;
	int t;
	int v;

	DgaBasisMon(const char* str_mon, bool known_diff_, PolyRaw&& diff_, int ss, int tt, int vv) :
		mon(str_to_mon(str_mon)), diff_determined(known_diff_), diff(diff_), s(ss), t(tt), v(vv) {};
	DgaBasisMon(std::vector<int> mon_, bool known_diff_, PolyRaw&& diff_, int ss, int tt, int vv) :
		mon(mon_), diff_determined(known_diff_), diff(diff_), s(ss), t(tt), v(vv) {};
};

PolyRaw str_to_poly(const char* str_poly)
{
	PolyRaw result;
	if (str_poly[0] == '\0')
		return std::vector<std::vector<int>>();
	std::stringstream ss(str_poly);
	while (ss.good()) {
		int i;
		ss >> i;
		if (result.empty())
			result.push_back(std::vector<int>());
		result.back().push_back(i);
		if (ss.peek() == ',')
			ss.ignore();
		else if (ss.peek() == ';') {
			ss.ignore();
			result.push_back(std::vector<int>());
		}
	}
	return result;
}

void load_generators(sqlite3* conn, const char* table_name, std::vector<DgaGenerator>& generators)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT gen_id, gen_name, gen_diff, s, t, v FROM ") + table_name + ";";
	sqlite3_prepare_v2(conn, cmd.c_str(), cmd.length() + 1, &stmt, NULL);
	while (sqlite3_step(stmt) == SQLITE_ROW)
		generators.emplace_back(sqlite3_column_int(stmt, 0), sqlite3_column_str(stmt, 1), str_to_poly(sqlite3_column_str(stmt, 2)),
			sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4), sqlite3_column_int(stmt, 5));
	sqlite3_finalize(stmt);
}

void load_relations(sqlite3* conn, const char* table_name, PolysRaw& relations)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT leading_term, basis FROM ") + table_name + ";";
	sqlite3_prepare_v2(conn, cmd.c_str(), cmd.length() + 1, &stmt, NULL);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		std::string str_leading_term = sqlite3_column_str(stmt, 0);
		std::string str_basis = sqlite3_column_str(stmt, 1);
		relations.push_back(str_to_poly((!str_basis.empty() ? str_leading_term + ';' + str_basis : str_leading_term).c_str()));
	}
	sqlite3_finalize(stmt);
}

void load_basis(sqlite3* conn, const char* table_name, std::vector<DgaBasisMon>& basis)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT mon, diff, s, t, v FROM ") + table_name + " WHERE t < 30;";
	sqlite3_prepare_v2(conn, cmd.c_str(), cmd.length() + 1, &stmt, NULL);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		if (sqlite3_column_type(stmt, 1) == SQLITE_NULL)
			basis.emplace_back(sqlite3_column_str(stmt, 0), false, PolyRaw(),
				sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4));
		else
			basis.emplace_back(sqlite3_column_str(stmt, 0), true, str_to_poly(sqlite3_column_str(stmt, 1)),
				sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4));
	}
	sqlite3_finalize(stmt);
}

/* num_mons[t] is the number of monomials up to degree t */
void load_mon_indices_tsv(sqlite3* conn, const char* table_name, std::vector<std::array<int, 5>>& num_mons_tsv)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT t, s, v, a, b FROM ") + table_name + " ORDER BY t, s, v;";
	sqlite3_prepare_v2(conn, cmd.c_str(), cmd.length() + 1, &stmt, NULL);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		num_mons_tsv.push_back({ sqlite3_column_int(stmt, 0), sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4) });
	}
	sqlite3_finalize(stmt);
}

inline std::pair<int, int> get_ab(const std::vector<std::array<int, 5>>& indices_tsv, int s, int t, int v)
{
	auto row = std::lower_bound(indices_tsv.begin(), indices_tsv.end(), std::array<int, 5>({ t, s, v, 0, 0 }));
	return std::pair<int, int>({ (*row)[3], (*row)[4] });
}

inline size_t get_index(const std::vector<DgaBasisMon>& basis, const MonRaw& mon, int a, int b)
{
	auto index = std::lower_bound(basis.begin() + a, basis.begin() + b, mon, [](const DgaBasisMon& b, const MonRaw& mon) { return cmp_mons(b.mon, mon); });
	return index - basis.begin();
}

void generate_diff(sqlite3* conn, const char* table_name_basis, const std::vector<DgaGenerator>& generators,
	const PolysRaw& relations)
{
	int v_diff = -2;
	/* Compile SQL statements */
	sqlite3_stmt* stmt_update_basis;
	std::string cmd_update_basis = std::string("UPDATE ") + table_name_basis + " SET diff=?1 WHERE mon_id=?2;";
	sqlite3_prepare_v2(conn, cmd_update_basis.c_str(), cmd_update_basis.length() + 1, &stmt_update_basis, NULL);

	/* Load the basis and indices. num_mons_t[t] is the number of monomials in degree <= t.
	In num_mons_tsv[i] = (t, s, v, a, b), [a, b) is the range of indices of monomial with degree (t, s, v). */
	std::vector<int> num_mons_t;
	std::string table_indices_t = std::string(table_name_basis) + "_indices_t";
	load_mon_indices(conn, table_indices_t.c_str(), num_mons_t);
	std::cout << "num_mon_t loaded! Size=" << num_mons_t.size() << '\n';
	std::vector<std::array<int, 5>> indices_tsv;
	std::string table_indices_tsv = std::string(table_name_basis) + "_indices_stv";
	load_mon_indices_tsv(conn, table_indices_tsv.c_str(), indices_tsv);
	std::cout << "num_mons_tsv loaded! Size=" << indices_tsv.size() << '\n';
	std::vector<DgaBasisMon> basis;
	load_basis(conn, table_name_basis, basis);
	std::cout << "basis loaded!\n";

	int prev_t = 0;
	execute_cmd(conn, "BEGIN TRANSACTION");
	for (size_t i = 0; i < basis.size(); i++) {
		if (!basis[i].diff_determined) {
			if (basis[i].t > prev_t) {
				prev_t = basis[i].t;
				std::cout << "t=" << basis[i].t << '\n';
				execute_cmd(conn, "END TRANSACTION");
				execute_cmd(conn, "BEGIN TRANSACTION");
			}
			if (basis[i].mon.empty()) {
				basis[i].diff = PolyRaw();
				basis[i].diff_determined = true;
			}
			else {
				int gen_id = basis[i].mon.front();
				MonRaw mon1 = div(basis[i].mon, { gen_id, 1 });
				auto ab = get_ab(indices_tsv, basis[i].s - generators[gen_id].s, basis[i].t - generators[gen_id].t, basis[i].v - generators[gen_id].v);
				size_t index_mon1 = get_index(basis, mon1, ab.first, ab.second);
				basis[i].diff = reduce(add(mul(generators[gen_id].diff, mon1), mul({ { gen_id, 1 } }, basis[index_mon1].diff)), relations);
				for (const auto& m : basis[i].diff) {
					auto ab_diff = get_ab(indices_tsv, basis[i].s + 1, basis[i].t, basis[i].v - 2);
					size_t index_diff = get_index(basis, m, ab_diff.first, ab_diff.second);
					basis[i].diff_indices.push_back(index_diff - ab_diff.first);
				}
				basis[i].diff_determined = true;
			}
			std::string str_diff(mon_to_str(basis[i].diff_indices));
			std::cout << "d(" << basis[i].mon << ")=" << str_diff << '\n';
			sqlite3_bind_text(stmt_update_basis, 1, str_diff.c_str(), str_diff.size() + 1, NULL);
			sqlite3_bind_int(stmt_update_basis, 2, i);
			sqlite3_step(stmt_update_basis);
			sqlite3_reset(stmt_update_basis);
		}
	}
	execute_cmd(conn, "END TRANSACTION");
	sqlite3_finalize(stmt_update_basis);
}

int main_test()
{
	PolyRaw rel1 = { {1, 1, 2, 1}, {0, 1} };
	PolyRaw rel2 = { {3, 1}, {0, 1} };
	PolysRaw rels = { rel1, rel2 };

	PolyRaw p = { {2, 1, 3, 1}, {0, 1, 1, 1, 2, 1} };
	std::cout << reduce(p, rels) << '\n';
	return 0;
}

int main_generate_diff(int argc, char** argv)
{
	//return main_test();

	sqlite3* conn;
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\ss.db)", &conn);

	const char* table_generators, * table_relations, * table_basis;
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
	std::vector<DgaGenerator> generators;
	load_generators(conn, table_generators, generators);

	PolysRaw relations;
	load_relations(conn, table_relations, relations);

	generate_diff(conn, table_basis, generators, relations);

	sqlite3_close(conn);
	return 0;
}