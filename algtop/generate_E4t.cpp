#include "database.h"

void load_reprs(sqlite3* conn, const std::string& table_name, array3d& reprs)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT repr FROM ") + table_name + ";";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		reprs.push_back(str_to_array2d(sqlite3_column_str(stmt, 0)));
	}
	sqlite3_finalize(stmt);
}

/* m1 is a sparse vector while [p1, p2) is a dense vector starting with index `index` */
bool divides(const array& m, array::const_iterator p1, array::const_iterator p2, size_t index)
{
	for (size_t k = 0; k < m.size(); k += 2) {
		if (m[k + 1] > *(p1 + (size_t(m[k]) - index)))
			return false;
	}
	return true;
}

array2d basis(const array3d& leadings, const std::vector<Deg>& gen_degs, Deg deg)
{
	array2d result;
	array mon_dense;
	mon_dense.resize(gen_degs.size());
	size_t index = mon_dense.size() - 1;

	while (true) {
		//std::cout << "index=" << index << '\n';
		mon_dense[index] = T_MAX;
		int e_max = T_MAX;
		for (const array& lead : leadings[index]) {
			if (divides(lead, mon_dense.begin() + index, mon_dense.end(), index)) {
				if (lead[1] - 1 < e_max)
					e_max = lead[1] - 1;
			}
		}
		int qs = gen_degs[index].s ? deg.s / gen_degs[index].s : T_MAX;
		int qt = gen_degs[index].t ? deg.t / gen_degs[index].t : T_MAX;
		int qv = gen_degs[index].v ? deg.v / gen_degs[index].v : T_MAX;
		e_max = std::min({ e_max, qs, qt, qv });
		//std::cout << "e_max=" << e_max << '\n';
		mon_dense[index] = e_max;
		deg -= gen_degs[index] * e_max;

		bool move_right = false;
		if (deg != Deg{ 0, 0, 0 }) {
			if (index > 0)
				index--;
			else
				move_right = true;
		}
		else {
			array mon_sparse;
			for (size_t i = index; i < mon_dense.size(); ++i) {
				if (mon_dense[i]) {
					mon_sparse.push_back(int(i));
					mon_sparse.push_back(mon_dense[i]);
				}
			}
			result.push_back(std::move(mon_sparse));
			if (index > 0) {
				deg += gen_degs[index];
				mon_dense[index--]--;
			}
			else
				move_right = true;
		}
		if (move_right) {
			for (deg += gen_degs[index] * mon_dense[index], ++index; index < mon_dense.size() && mon_dense[index] == 0; ++index);
			if (index == mon_dense.size()) {
				//std::cout << "mon_dense=" << mon_dense << '\n';
				break;
			}
			else {
				deg += gen_degs[index];
				mon_dense[index--]--;
			}
		}
	}
	return result;
}

/* Assume poly is a boundary. Return the chain with it as boundary */
array2d d_inv(const array2d& poly, const array3d& gb, const array3d& leadings, const std::vector<Deg>& gen_degs, const array3d& diffs)
{
	if (poly.empty())
		return array2d();
	Deg deg_poly = deg(poly[0], gen_degs);
	Deg deg_result = deg_poly - Deg{ 1, 0, -2 };
	array2d basis_in_poly = basis(leadings, gen_degs, deg_poly);
	array2d basis_in_result = basis(leadings, gen_degs, deg_result);
	std::sort(basis_in_poly.begin(), basis_in_poly.end(), cmp_mons);
	std::sort(basis_in_result.begin(), basis_in_result.end(), cmp_mons);
	array2d map_d;
	for (const array& mon : basis_in_result)
		map_d.push_back(poly_to_indices(reduce(get_diff(mon, diffs), gb), basis_in_poly));
	array2d image, kernel, g;
	set_linear_map(map_d, image, kernel, g);
	return indices_to_poly(get_image(image, g, poly_to_indices(poly, basis_in_poly)), basis_in_result);
}

void generate_E4bk(sqlite3* conn, const std::string& table_prefix, const std::string& table1_prefix, const array2d& a, int bk_gen_id, int t_max)
{
	/* Load generators */

	std::vector<Deg> gen_degs;
	load_gen_degs(conn, table_prefix + "_generators", gen_degs);
	std::cout << "gen_degs.size()=" << gen_degs.size() << '\n';

	array3d reprs;
	load_reprs(conn, table_prefix + "_generators", reprs);
	std::cout << "reprs.size()=" << reprs.size() << '\n';

	std::vector<Deg> gen_degs_E2t;
	load_gen_degs(conn, "E2t_generators", gen_degs_E2t);
	while (gen_degs_E2t.size() > size_t(bk_gen_id) + 1) gen_degs_E2t.pop_back();
	std::cout << "gen_degs_E2t.size()=" << gen_degs_E2t.size() << '\n';

	array3d diffs_E2t;
	load_gen_diffs(conn, "E2t_generators", diffs_E2t);
	std::cout << "diffs_E2t.size()=" << diffs_E2t.size() << '\n';

	array gen_degs_t;
	for (auto p = gen_degs.begin(); p < gen_degs.end(); ++p)
		gen_degs_t.push_back(p->t);

	/* Load Groebner basis */

	array3d gb;
	load_gb(conn, table_prefix + "_relations", gb);
	std::cout << "gb.size()=" << gb.size() << '\n';

	array3d gb_E2t;
	load_gb(conn, "E2_relations", gb_E2t);
	std::cout << "gb_E2t.size()=" << gb_E2t.size() << '\n';

	array3d leadings; // TODO: save greobner basis as array4d
	leadings.resize(gen_degs.size());
	for (const array2d& g : gb)
		leadings[g[0][0]].push_back(g[0]);

	array3d leadings_E2t;
	leadings_E2t.resize(gen_degs_E2t.size());
	for (const array2d& g : gb_E2t)
		leadings_E2t[g[0][0]].push_back(g[0]);

	/* Compute new generators */

	array4d b_ = ann_seq(gb, { a }, gen_degs_t, t_max);
	std::cout << "b_=" << b_ << '\n';
	std::cout << "b_.size()=" << b_.size() << '\n';

	array3d b;
	for (array3d& v : b_)
		b.push_back(std::move(v[0]));
	array4d c_ = ann_seq(gb, b, gen_degs_t, t_max - gen_degs_E2t[bk_gen_id].t);
	std::cout << "c_.size()=" << c_.size() << '\n';

	std::vector<Deg> gen_degs_E2bk_1(gen_degs_E2t);
	gen_degs_E2bk_1.pop_back();
	array3d ab_inv;
	for (const array2d& bi : b) {
		array2d abi = mul(a, bi);
		array2d abi_repr = evaluate(abi, [&reprs](int i) {return reprs[i]; }, gb_E2t);
		ab_inv.push_back(d_inv(abi_repr, gb_E2t, leadings_E2t, gen_degs_E2bk_1, diffs_E2t));//
	}

	if (gen_degs_E2t.back().t * 2 <= t_max) { /* Add [x^2] */
		gen_degs.push_back(gen_degs_E2t.back() * 2);
		gen_degs_t.push_back(gen_degs.back().t);
		reprs.push_back({ {bk_gen_id, 2} });
	}
	for (size_t i = 0; i < b.size(); ++i) { /* Add [xb+d^{-1}(ab)] */
		gen_degs.push_back(gen_degs_E2t[bk_gen_id] + deg(b[i], gen_degs));
		gen_degs_t.push_back(gen_degs.back().t);
		reprs.push_back(add(mul(evaluate(b[i], [&reprs](int i) {return reprs[i]; }, gb_E2t), { bk_gen_id, 1 }), ab_inv[i]));
	}

	/* Compute the relations */

	add_rels(gb, { a }, gen_degs_t, t_max); /* Add relations dx=0 */
	leadings.clear();
	leadings.resize(gen_degs.size());
	for (const array2d& g : gb)
		leadings[g[0][0]].push_back(g[0]);

	std::vector<Deg> rel_degs;
	for (size_t i = gen_degs.size() - b.size(); i < gen_degs.size(); ++i)
		for (size_t j = i; j < gen_degs.size(); ++j) {
			Deg deg_gigj = gen_degs[i] + gen_degs[j];
			if (deg_gigj.t <= t_max)
				rel_degs.push_back(deg_gigj);
		}
	for (const array3d& ci : c_)
		for (size_t j = 0; j < b.size(); ++j)
			if (!ci[j].empty()) {
				rel_degs.push_back(gen_degs[gen_degs.size() - b.size() + j] + deg(ci[j], gen_degs));
				break;
			}
	std::sort(rel_degs.begin(), rel_degs.end());

	array3d rels;
	for (size_t i = 0; i < rel_degs.size(); ++i) {
		const Deg d = rel_degs[i];
		std::cout << i << '/' << rel_degs.size() << ' ' << "deg=" << '(' << d.s << ", " << d.t << ", " << d.v << ')' << "          \r";
		array2d basis_d = basis(leadings, gen_degs, d);
		array2d basis_d_E2t = basis(leadings_E2t, gen_degs_E2t, d);
		array2d basis_d_E2t1 = basis(leadings_E2t, gen_degs_E2t, d - Deg{ 1, 0, -2 });
		std::sort(basis_d.begin(), basis_d.end(), cmp_mons);
		std::sort(basis_d_E2t.begin(), basis_d_E2t.end(), cmp_mons);
		std::sort(basis_d_E2t1.begin(), basis_d_E2t1.end(), cmp_mons);

		array2d map_diff;
		for (const array& mon : basis_d_E2t1)
			map_diff.push_back(poly_to_indices(reduce(get_diff(mon, diffs_E2t), gb_E2t), basis_d_E2t));
		array2d image_diff, kernel_diff, g_diff;
		set_linear_map(map_diff, image_diff, kernel_diff, g_diff);

		array2d map_repr;
		for (const array& mon : basis_d) {
			array repr = residue(image_diff, poly_to_indices(evaluate({ mon }, [&reprs](int i) {return reprs[i]; }, gb_E2t), basis_d_E2t));
			map_repr.push_back(std::move(repr));
		}
		array2d image_repr, kernel_repr, g_repr;
		set_linear_map(map_repr, image_repr, kernel_repr, g_repr);
		for (const array& rel_indices : kernel_repr)
			rels.push_back(indices_to_poly(rel_indices, basis_d));
	}

	size_t old_gb_size = gb.size();
	add_rels(gb, rels, gen_degs_t, t_max);
	std::cout << "Groebner size increases by " << gb.size() - old_gb_size << '\n';


	/* Save generators */
	sqlite3_stmt* stmt_update_generators;
	std::string cmd_update_generators = std::string("INSERT INTO ") + table1_prefix + "_generators (gen_id, repr, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5);";
	sqlite3_prepare_v100(conn, cmd_update_generators, &stmt_update_generators);

	execute_cmd(conn, "BEGIN TRANSACTION");
	for (size_t i = 0; i < gen_degs.size(); ++i) {
		sqlite3_bind_int(stmt_update_generators, 1, int(i));
		sqlite3_bind_str(stmt_update_generators, 2, array2d_to_str(reprs[i]));
		sqlite3_bind_int(stmt_update_generators, 3, gen_degs[i].s);
		sqlite3_bind_int(stmt_update_generators, 4, gen_degs[i].t);
		sqlite3_bind_int(stmt_update_generators, 5, gen_degs[i].v);
		sqlite3_step(stmt_update_generators);
		sqlite3_reset(stmt_update_generators);
	}
	execute_cmd(conn, "END TRANSACTION");
	std::cout << gen_degs.size() << " Generators are inserted!\n";

	sqlite3_finalize(stmt_update_generators);

	/* Save relations */
	sqlite3_stmt* stmt_update_relations;
	std::string cmd_update_relations = std::string("INSERT INTO ") + table1_prefix + "_relations (leading_term, basis, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5);"; // Insert s, t, v
	sqlite3_prepare_v100(conn, cmd_update_relations, &stmt_update_relations);

	execute_cmd(conn, "BEGIN TRANSACTION");
	for (size_t i = 0; i < gb.size(); ++i) {
		Deg d = deg(gb[i], gen_degs);
		sqlite3_bind_str(stmt_update_relations, 1, array_to_str(gb[i][0]));
		sqlite3_bind_str(stmt_update_relations, 2, array2d_to_str(gb[i].begin() + 1, gb[i].end()));
		sqlite3_bind_int(stmt_update_relations, 3, d.s);
		sqlite3_bind_int(stmt_update_relations, 4, d.t);
		sqlite3_bind_int(stmt_update_relations, 5, d.v);
		sqlite3_step(stmt_update_relations);
		sqlite3_reset(stmt_update_relations);
	}
	execute_cmd(conn, "END TRANSACTION");
	std::cout << gb.size() << " relations are inserted!\n";

	sqlite3_finalize(stmt_update_relations);
}

int main_test1(int argc, char** argv)
{
	//array a = { 1, 2, 3, 4 };
	//std::cout << std::hash<array>(a);
	sqlite3* conn;
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\ss.db)", &conn);

	/* Load basis */

	std::map<Deg, array2d> basis;
	load_basis(conn, "E2_basis", basis);
	for (const auto& pair : basis) {
		std::cout << pair.first << ": " << pair.second.size() << '\n';
	}

	sqlite3_close(conn);
	return 0;
}

int main_generate_E4t(int argc, char** argv)
{
	//return main_test1(argc, argv);

	sqlite3* conn;
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)", &conn);

	int t_max = 74;
	std::cout << "E4b1\n";
	generate_E4bk(conn, "E4", "E4b1", { {0, 1} }, 58, t_max);
	std::cout << "E4b2\n";
	generate_E4bk(conn, "E4b1", "E4b2", { {54, 1} }, 59, t_max);
	std::cout << "E4b3\n";
	generate_E4bk(conn, "E4b2", "E4b3", { {75, 1} }, 60, t_max);
	std::cout << "E4b4\n";
	generate_E4bk(conn, "E4b3", "E4b4", { {113, 1} }, 61, t_max);
	std::cout << "E4b5\n";
	generate_E4bk(conn, "E4b4", "E4b5", { {172, 1} }, 62, t_max);
	std::cout << "E4b6\n";
	generate_E4bk(conn, "E4b5", "E4b6", { {249, 1} }, 63, t_max);
	//std::cout << "E4b7\n";
	//generate_E4bk(conn, "E4b6", "E4b7", 64);

	sqlite3_close(conn);
	return 0;
}