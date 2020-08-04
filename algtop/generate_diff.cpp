#include "database.h"
#include <sstream>


struct DgaGenerator
{
	int id;
	std::string name;
	array2d diff;
	int s;
	int t;
	int v;

	DgaGenerator(int gen_id, const char* gen_name, array2d&& gen_diff, int s_, int t_, int v_) :
		id(gen_id), name(gen_name), diff(gen_diff), s(s_), t(t_), v(v_) {};
};

struct DgaBasisMon1
{
	array mon;
	bool diff_determined;
	array2d diff;
	array diff_indices;
	int s;
	int t;
	int v;

	DgaBasisMon1(const char* str_mon, bool diff_determined_, array&& diff_indices_, int s_, int t_, int v_) :
		mon(str_to_array(str_mon)), diff_determined(diff_determined_), diff_indices(diff_indices_), s(s_), t(t_), v(v_) {};
	DgaBasisMon1(array&& mon_, bool diff_determined_, array&& diff_indices_, int s_, int t_, int v_) :
		mon(mon_), diff_determined(diff_determined_), diff_indices(diff_indices_), s(s_), t(t_), v(v_) {};
};

inline array2d indices_to_poly(const array& indices, const std::vector<DgaBasisMon1>& basis, int a)
{
	array2d result;
	for (int i : indices)
		result.push_back(basis[size_t(a) + i].mon);
	return result;
}

inline int get_index(const std::vector<DgaBasisMon1>& basis, const array& mon, int a, int b)
{
	auto index = std::lower_bound(basis.begin() + a, basis.begin() + b, mon, [](const DgaBasisMon1& b, const array& mon) { return cmp_mons(b.mon, mon); });
	return int(index - basis.begin());
}

inline array poly_to_indices(const array2d& poly, const std::vector<DgaBasisMon1>& basis, int a, int b)
{
	array result;
	for (const array& mon : poly)
		result.push_back(get_index(basis, mon, a, b));
	return result;
}

void load_generators(sqlite3* conn, const char* table_name, std::vector<DgaGenerator>& generators)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT gen_id, gen_name, gen_diff, s, t, v FROM ") + table_name + ";";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW)
		generators.emplace_back(sqlite3_column_int(stmt, 0), sqlite3_column_str(stmt, 1), str_to_array2d(sqlite3_column_str(stmt, 2)),
			sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4), sqlite3_column_int(stmt, 5));
	sqlite3_finalize(stmt);
}

void load_relations(sqlite3* conn, const std::string& table_name, array3d& relations)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT leading_term, basis FROM ") + table_name + ";";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		std::string str_leading_term = sqlite3_column_str(stmt, 0);
		std::string str_basis = sqlite3_column_str(stmt, 1);
		relations.push_back(str_to_array2d((!str_basis.empty() ? str_leading_term + ';' + str_basis : str_leading_term).c_str()));
	}
	sqlite3_finalize(stmt);
}

void load_basis(sqlite3* conn, const char* table_name, std::vector<DgaBasisMon1>& basis, const std::vector<std::array<int, 5>>& indices_tsv, int r)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT mon, diff, s, t, v FROM ") + table_name + " ;";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	int prev_t = 0;
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		int t = sqlite3_column_int(stmt, 3);
		if (t > prev_t) {
			std::cout << "t=" << t << '\n';
			prev_t = t;
		}
		if (sqlite3_column_type(stmt, 1) == SQLITE_NULL) {
			basis.emplace_back(sqlite3_column_str(stmt, 0), false, array(),
				sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4));
		}
		else {
			basis.emplace_back(sqlite3_column_str(stmt, 0), true, str_to_array(sqlite3_column_str(stmt, 1)),
				sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4));
		}
	}
	for (auto& base : basis) {
		auto ab = get_ab(indices_tsv, base.s + 1, base.t, base.v - r);
		base.diff = indices_to_poly(base.diff_indices, basis, ab.first);
	}
	sqlite3_finalize(stmt);
}

/* num_mons[t] is the number of monomials up to degree t */
void load_mon_indices_tsv(sqlite3* conn, const std::string& table_name, std::vector<std::array<int, 5>>& num_mons_tsv)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT t, s, v, a, b FROM ") + table_name + " ORDER BY t, s, v;";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		num_mons_tsv.push_back({ sqlite3_column_int(stmt, 0), sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4) });
	}
	sqlite3_finalize(stmt);
}

void generate_diff(sqlite3* conn, const char* table_name_basis, const std::vector<DgaGenerator>& generators,
	const array3d& relations, int r)
{
	/* Compile SQL statements */
	sqlite3_stmt* stmt_update_basis;
	std::string cmd_update_basis = std::string("UPDATE ") + table_name_basis + " SET diff=?1 WHERE mon_id=?2;";
	sqlite3_prepare_v100(conn, cmd_update_basis, &stmt_update_basis);

	/* Load the basis and indices. num_mons_t[t] is the number of monomials in degree <= t.
	In num_mons_tsv[i] = (t, s, v, a, b), [a, b) is the range of indices of monomial with degree (t, s, v). */
	array num_mons_t;
	std::string table_indices_t = std::string(table_name_basis) + "_indices_t";
	load_mon_indices(conn, table_indices_t.c_str(), num_mons_t);
	std::cout << "num_mon_t loaded! Size=" << num_mons_t.size() << '\n';

	std::vector<std::array<int, 5>> indices_tsv;
	std::string table_indices_tsv = std::string(table_name_basis) + "_indices_stv";
	load_mon_indices_tsv(conn, table_indices_tsv.c_str(), indices_tsv);
	std::cout << "num_mons_tsv loaded! Size=" << indices_tsv.size() << '\n';

	std::vector<DgaBasisMon1> basis;
	load_basis(conn, table_name_basis, basis, indices_tsv, r);
	std::cout << "basis loaded!\n";

	execute_cmd(conn, "BEGIN TRANSACTION");
	int prev_t = 0;
	for (size_t i = 0; i < basis.size(); i++) {
		if (!basis[i].diff_determined) {
			array diff_indices;
			if (basis[i].t > prev_t) {
				prev_t = basis[i].t;
				std::cout << "t=" << basis[i].t << '\n';
				execute_cmd(conn, "END TRANSACTION");
				execute_cmd(conn, "BEGIN TRANSACTION");
			}
			if (basis[i].mon.empty()) {
				basis[i].diff = array2d();
				basis[i].diff_determined = true;
			}
			else {
				int gen_id = basis[i].mon.front();
				array mon1 = div(basis[i].mon, { gen_id, 1 });
				auto ab = get_ab(indices_tsv, basis[i].s - generators[gen_id].s, basis[i].t - generators[gen_id].t, basis[i].v - generators[gen_id].v);
				size_t index_mon1 = get_index(basis, mon1, ab.first, ab.second);
				basis[i].diff = reduce(add(mul(generators[gen_id].diff, mon1), mul({ { gen_id, 1 } }, basis[index_mon1].diff)), relations);
				auto ab_diff = get_ab(indices_tsv, basis[i].s + 1, basis[i].t, basis[i].v - r);
				diff_indices = poly_to_indices(basis[i].diff, basis, ab_diff.first, ab_diff.second);
				basis[i].diff_determined = true;
			}
			std::string str_diff(array_to_str(diff_indices));
			//std::cout << "d(" << basis[i].mon << ")=" << str_diff << '\n';
			sqlite3_bind_str(stmt_update_basis, 1, str_diff);
			sqlite3_bind_int(stmt_update_basis, 2, int(i));
			sqlite3_step(stmt_update_basis);
			sqlite3_reset(stmt_update_basis);
		}
	}
	execute_cmd(conn, "END TRANSACTION");
	sqlite3_finalize(stmt_update_basis);
}

int main_test1()
{
	std::vector<DgaBasisMon1> basis;
	std::cout << basis.max_size() << '\n';
	for (int i = 0; i < 20000000; i++) {
		if (i % 1000000 == 0)
			std::cout << i << '\n';
		basis.emplace_back(array({ 0, 38, 1, 2 }), false, array(),
			0, 0, 0);
	}
	return 0;
}

int main_generate_diff(int argc, char** argv)
{
	//return main_test1();

	sqlite3* conn;
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\ss.db)", &conn);

	const char* table_generators, * table_relations, * table_basis;
	if (argc == 1) {
		table_generators = "E2_generators";
		table_relations = "E2_relations";
		table_basis = "E2_basis";
	}
	else {
		table_generators = argv[1];
		table_relations = argv[2];
		table_basis = argv[3];
	}
	std::vector<DgaGenerator> generators;
	load_generators(conn, table_generators, generators);

	array3d relations;
	load_relations(conn, table_relations, relations);

	generate_diff(conn, table_basis, generators, relations, 2);

	sqlite3_close(conn);
	return 0;
}