#include "main.h"
#include <chrono>
#include <fstream>
#include <future>

/********** STRUCTS AND CLASSES **********/

struct DgaBasis1
{
	Mon1d basis;
	Poly1d diffs;
};

/********** FUNCTIONS **********/

void load_dga_basis(const Database& db, const std::string& table_name, std::map<Deg, DgaBasis1>& basis, int r)
{
	Statement stmt;
	stmt.init(db, "SELECT mon, diff, s, t, v FROM " + table_name + " ORDER BY mon_id;");
	int prev_t = 0;
	std::map<Deg, array2d> mon_diffs;
	while (stmt.step() == SQLITE_ROW) {
		Deg d = { stmt.column_int(2), stmt.column_int(3), stmt.column_int(4) };
		if (d.t > prev_t) {
			std::cout << "load_dag_basis, t=" << d.t << "          \r";
			prev_t = d.t;
		}
		basis[d].basis.push_back(str_to_Mon(stmt.column_str(0)));
		basis[d].diffs.emplace_back();
		mon_diffs[d].push_back({ str_to_array(stmt.column_str(1)) });
	}
	for (auto& [d, basis_d] : basis) {
		for (size_t i = 0; i < basis_d.basis.size(); ++i) {
			basis_d.diffs[i] = indices_to_Poly(mon_diffs[d][i], basis[d + Deg{ 1, 0, -r }].basis);
		}
	}
}

/* m1 is a sparse vector while [p1, p2) is a dense vector starting with index `index` */
bool divides(const Mon& m, array::const_iterator p1, array::const_iterator p2, size_t index)
{
	for (MonInd p = m.begin(); p != m.end(); ++p)
		if (p->exp > *(p1 + (size_t(p->gen) - index)))
			return false;
	return true;
}

/* return a basis in degree `deg` */
Mon1d get_basis(const Mon2d& leadings, const std::vector<Deg>& gen_degs, Deg deg)
{
	Mon1d result;
	array mon_dense;
	mon_dense.resize(gen_degs.size());
	size_t index = mon_dense.size() - 1;

	while (true) {
		//std::cout << "index=" << index << '\n';
		mon_dense[index] = T_MAX;
		int e_max = T_MAX;
		for (const Mon& lead : leadings[index]) {
			if (divides(lead, mon_dense.begin() + index, mon_dense.end(), index)) {
				if (lead[0].exp - 1 < e_max)
					e_max = lead[0].exp - 1;
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
			Mon mon;
			for (size_t i = index; i < mon_dense.size(); ++i)
				if (mon_dense[i])
					mon.emplace_back(int(i), mon_dense[i]);
			result.push_back(std::move(mon));
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

/* return a basis for polynomails of bi for t<=t_max */
std::map<Deg, DgaBasis1> get_basis_bi(const std::vector<Deg>& gen_degs, const Poly1d& gen_diffs, int t_max)
{
	std::map<Deg, DgaBasis1> result;
	result[Deg{ 0, 0, 0 }].basis.push_back({});
	result[Deg{ 0, 0, 0 }].diffs.push_back({});

	for (int t = 1; t <= t_max; ++t) {
		std::map<Deg, DgaBasis1> basis_new;
		std::cout << "Computing basis_bi, t=" << t << "          \r";
		for (int gen_id = (int)gen_degs.size() - 1; gen_id >= 58; --gen_id) { /* 58 is the gen_id of b1 */
			int t1 = t - gen_degs[gen_id].t;
			if (t1 >= 0) {
				auto p1 = result.lower_bound(Deg{ 0, t1, 0 });
				auto p2 = result.lower_bound(Deg{ 0, t1 + 1, 0 });
				for (auto p = p1; p != p2; ++p) {
					auto p_m = p->second.basis.begin();
					auto p_diff = p->second.diffs.begin();
					for (; p_m != p->second.basis.end(); ++p_m, ++p_diff) {
						if (p_m->empty() || gen_id <= p_m->front().gen) {
							Mon mon(mul(*p_m, { {gen_id, 1} }));

							int s = p->first.s + gen_degs[gen_id].s;
							int v = p->first.v + gen_degs[gen_id].v;
							basis_new[Deg{ s, t, v }].basis.push_back(std::move(mon));
							Poly diff = add(mul(gen_diffs[gen_id], *p_m), mul({ {gen_id, 1} }, *p_diff));
							basis_new[Deg{ s, t, v }].diffs.push_back(std::move(diff));
						}
					}
				}
			}
		}
		result.merge(basis_new);
	}
	return result;
}

/* build the basis for t<=t_max */
void get_basis(const Mon2d& leadings, const std::vector<Deg>& gen_degs, std::map<Deg, Mon1d>& basis, int t_max)
{
	int t_min;
	if (basis.empty()) {
		basis[Deg{ 0, 0, 0 }].push_back({}); /* if no monomial present insert the unit */
		t_min = 1;
	}
	else
		t_min = basis.rbegin()->first.t + 1;

	for (int t = t_min; t <= t_max; ++t) {
		std::map<Deg, Mon1d> basis_new;
		for (int gen_id = (int)gen_degs.size() - 1; gen_id >= 0; --gen_id) { /* 58 is the gen_id of b1 */
			int t1 = t - gen_degs[gen_id].t;
			if (t1 >= 0) {
				auto p1 = basis.lower_bound(Deg{ 0, t1, 0 });
				auto p2 = basis.lower_bound(Deg{ 0, t1 + 1, 0 });
				for (auto p = p1; p != p2; ++p) {
					for (auto p_m = p->second.begin(); p_m != p->second.end(); ++p_m) {
						if (p_m->empty() || gen_id <= p_m->front().gen) {
							Mon mon(mul(*p_m, { {gen_id, 1} }));
							if ((size_t)gen_id >= leadings.size() || std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(),
								[&mon](const Mon& _m) { return divides(_m, mon); }))
								basis_new[p->first + gen_degs[gen_id]].push_back(std::move(mon));
						}
					}
				}
			}
		}
		basis.merge(basis_new);
	}
}

void get_basis_E2t(const std::map<Deg, DgaBasis1>& basis_E2, const std::map<Deg, DgaBasis1>& basis_bi, Mon1d* p_basis, Poly1d* p_diffs, const Deg& deg)
{
	for (auto p2 = basis_bi.rbegin(); p2 != basis_bi.rend(); ++p2) {
		const Deg& d2 = p2->first;
		if (d2.s <= deg.s && d2.t <= deg.t && d2.v <= deg.v) {
			const Deg d1 = deg - d2;
			auto p1 = basis_E2.find(d1);
			if (p1 != basis_E2.end()) {
				if (p_diffs == nullptr) {
					for (size_t i = 0; i < p1->second.basis.size(); ++i)
						for (size_t j = 0; j < p2->second.basis.size(); ++j)
							p_basis->push_back(mul(p1->second.basis[i], p2->second.basis[j]));
				}
				else
					for (size_t i = 0; i < p1->second.basis.size(); ++i)
						for (size_t j = 0; j < p2->second.basis.size(); ++j) {
							p_basis->push_back(mul(p1->second.basis[i], p2->second.basis[j]));
							Poly diff = add(mul(p1->second.diffs[i], p2->second.basis[j]), mul(p2->second.diffs[j], p1->second.basis[i]));
							p_diffs->push_back(std::move(diff));
						}
			}
		}
	}
}

/* Assume poly is a boundary. Return the chain with it as boundary */
Poly d_inv(const Poly& poly, const Poly1d& gb, const Mon2d& leadings, const std::vector<Deg>& gen_degs, const Poly1d& diffs)
{
	if (poly.empty())
		return {};
	Deg deg_poly = get_deg(poly[0], gen_degs);
	Deg deg_result = deg_poly - Deg{ 1, 0, -2 };
	Mon1d basis_in_poly = get_basis(leadings, gen_degs, deg_poly);
	Mon1d basis_in_result = get_basis(leadings, gen_degs, deg_result);
	std::sort(basis_in_poly.begin(), basis_in_poly.end());
	std::sort(basis_in_result.begin(), basis_in_result.end());
	array2d map_d;
	for (const Mon& mon : basis_in_result)
		map_d.push_back(Poly_to_indices(reduce(get_diff(mon, diffs), gb), basis_in_poly));
	array2d image, kernel, g;
	set_linear_map(map_d, image, kernel, g);
	return indices_to_Poly(get_image(image, g, Poly_to_indices(poly, basis_in_poly)), basis_in_result);
}

std::vector<rel_heap_t> find_relations(Deg d, Mon1d& basis_d, const std::map<Deg, DgaBasis1>& basis_E2, const std::map<Deg, DgaBasis1>& basis_bi, const Poly1d& gb_E2t, const Poly1d& reprs)
{
	Mon1d basis_d_E2t;
	get_basis_E2t(basis_E2, basis_bi, &basis_d_E2t, nullptr, d);
	Mon1d basis_d1_E2t;
	Poly1d diffs_d_E2t;
	get_basis_E2t(basis_E2, basis_bi, &basis_d1_E2t, &diffs_d_E2t, d - Deg{ 1, 0, -2 });
	array indices = range((int)basis_d1_E2t.size());
	std::sort(indices.begin(), indices.end(), [&basis_d1_E2t](int i1, int i2) {return basis_d1_E2t[i1] < basis_d1_E2t[i2]; });
	std::sort(basis_d_E2t.begin(), basis_d_E2t.end());
	std::sort(basis_d1_E2t.begin(), basis_d1_E2t.end());

	array2d map_diff;
	for (int i : indices)
		map_diff.push_back(Poly_to_indices(reduce(diffs_d_E2t[indices[i]], gb_E2t), basis_d_E2t));
	array2d image_diff, kernel_diff, g_diff;
	set_linear_map(map_diff, image_diff, kernel_diff, g_diff);

	array2d map_repr;
	for (const Mon& mon : basis_d) {
		array repr = residue(image_diff, Poly_to_indices(evaluate({ mon }, [&reprs](int i) {return reprs[i]; }, gb_E2t), basis_d_E2t));
		map_repr.push_back(std::move(repr));
	}
	array2d image_repr, kernel_repr, g_repr;
	set_linear_map(map_repr, image_repr, kernel_repr, g_repr);

	std::vector<rel_heap_t> result;
	for (const array& rel_indices : kernel_repr)
		result.push_back(rel_heap_t{ indices_to_Poly(rel_indices, basis_d), d.t });

	for (const array& rel_indices : kernel_repr)
		basis_d[rel_indices[0]].clear();
	basis_d.erase(std::remove_if(basis_d.begin(), basis_d.end(), [](const Mon& m) {return m.empty(); }), basis_d.end());

	return result;
}

void generate_E4bk(const Database& db, const std::string& table_prefix, const std::string& table1_prefix, const Poly& a, int bk_gen_id, int t_max)
{
	std::vector<Deg> gen_degs = db.load_gen_degs(table_prefix + "_generators");
	array gen_degs_t;
	for (auto p = gen_degs.begin(); p < gen_degs.end(); ++p)
		gen_degs_t.push_back(p->t);

	Poly1d reprs = db.load_gen_reprs(table_prefix + "_generators");

	std::vector<Deg> gen_degs_E2t = db.load_gen_degs("E2t_generators");
	while (gen_degs_E2t.size() > size_t(bk_gen_id) + 1)
		gen_degs_E2t.pop_back();

	Poly1d diffs_E2t = db.load_gen_diffs("E2t_generators");

	Poly1d gb = db.load_gb(table_prefix + "_relations");
	Mon2d leadings;
	leadings.resize(gen_degs.size());
	for (const Poly& g : gb)
		leadings[g[0][0].gen].push_back(g[0]);

	Poly1d gb_E2t = db.load_gb("E2_relations");
	Mon2d leadings_E2t;
	leadings_E2t.resize(gen_degs_E2t.size());
	for (const Poly& g : gb_E2t)
		leadings_E2t[g[0][0].gen].push_back(g[0]);

	/* load E2_basis */

	std::map<Deg, DgaBasis1> basis_E2;
	load_dga_basis(db, "E2_basis", basis_E2, 2);
	std::cout << "basis_E2 loaded! Size=" << basis_E2.size() << '\n';

	/* generate basis of polynomials of b_1,...,b_k */

	std::map<Deg, DgaBasis1> basis_bi = get_basis_bi(gen_degs_E2t, diffs_E2t, t_max);
	std::cout << "basis_bi loaded! Size=" << basis_bi.size() << '\n';

	/* compute new generators */

	Poly2d b_ = ann_seq(gb, { a }, gen_degs_t, t_max);
	std::cout << "b_=" << b_ << '\n';
	std::cout << "b_.size()=" << b_.size() << '\n';

	Poly1d b;
	for (Poly1d& v : b_)
		b.push_back(std::move(v[0]));
	Poly2d c_ = ann_seq(gb, b, gen_degs_t, t_max - gen_degs_E2t[bk_gen_id].t);
	std::cout << "c_.size()=" << c_.size() << '\n';

	std::vector<Deg> gen_degs_E2bk_1(gen_degs_E2t);
	gen_degs_E2bk_1.pop_back();
	Poly1d ab_inv;
	for (const Poly& bi : b) {
		Poly abi = mul(a, bi);
		Poly abi_repr = evaluate(abi, [&reprs](int i) {return reprs[i]; }, gb_E2t);
		ab_inv.push_back(d_inv(abi_repr, gb_E2t, leadings_E2t, gen_degs_E2bk_1, diffs_E2t));//
	}

	if (gen_degs_E2t.back().t * 2 <= t_max) { /* Add [x^2] */
		gen_degs.push_back(gen_degs_E2t.back() * 2);
		gen_degs_t.push_back(gen_degs.back().t);
		reprs.push_back({ {{bk_gen_id, 2}} });
	}
	for (size_t i = 0; i < b.size(); ++i) { /* Add [xb+d^{-1}(ab)] */
		gen_degs.push_back(gen_degs_E2t[bk_gen_id] + get_deg(b[i], gen_degs));
		gen_degs_t.push_back(gen_degs.back().t);
		reprs.push_back(add(mul(evaluate(b[i], [&reprs](int i) {return reprs[i]; }, gb_E2t), { {bk_gen_id, 1} }), ab_inv[i]));
	}

	/* compute the relations */

	add_rels(gb, { a }, gen_degs_t, t_max); /* Add relations dx=0 */
	leadings.clear();
	leadings.resize(gen_degs.size());
	for (const Poly& g : gb)
		leadings[g[0][0].gen].push_back(g[0]);

	std::vector<Deg> rel_degs;
	for (size_t i = gen_degs.size() - b.size(); i < gen_degs.size(); ++i)
		for (size_t j = i; j < gen_degs.size(); ++j) {
			Deg deg_gigj = gen_degs[i] + gen_degs[j];
			if (deg_gigj.t <= t_max)
				rel_degs.push_back(deg_gigj);
		}
	for (const Poly1d& ci : c_)
		for (size_t j = 0; j < b.size(); ++j)
			if (!ci[j].empty()) {
				rel_degs.push_back(gen_degs[gen_degs.size() - b.size() + j] + get_deg(ci[j], gen_degs));
				break;
			}
	std::sort(rel_degs.begin(), rel_degs.end());
	rel_degs.erase(std::unique(rel_degs.begin(), rel_degs.end()), rel_degs.end());

	std::ofstream myfile;
	myfile.open("rel_gens.txt");
	for (const auto& d : rel_degs)
		myfile << d << '\n';
	myfile.close();

	std::map<Deg, Mon1d> basis;
	size_t i = 0;
	std::vector<rel_heap_t> heap;
	for (int t = 1; t <= t_max; ++t) {
		get_basis(leadings, gen_degs, basis, t);

		std::vector<std::future<std::vector<rel_heap_t>>> futures;
		for (; i < rel_degs.size() && rel_degs[i].t == t; ++i) {
			const Deg d = rel_degs[i];
			Mon1d& basis_d = basis[d];
			futures.push_back(std::async(std::launch::async, find_relations, d, std::ref(basis_d), std::ref(basis_E2), std::ref(basis_bi), std::ref(gb_E2t), std::ref(reprs)));
		}
		for (size_t j = 0; j < futures.size(); ++j) {
			futures[j].wait();
			std::cout << "t=" << t << " completed thread=" << j + 1 << '/' << futures.size() << "          \r";
		}
		for (auto& f : futures)
			for (auto& rel : f.get()) {
				heap.push_back(std::move(rel));
				std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
			}
		
		add_rels(gb, heap, gen_degs_t, std::min(t + 1, t_max), t_max);
		leadings.clear();
		leadings.resize(gen_degs.size());
		for (const Poly& g : gb)
			leadings[g[0][0].gen].push_back(g[0]);
	}

	/* Save generators */
	db.save_generators(table1_prefix + "_generators", gen_degs, reprs);

	/* Save relations */
	db.save_gb(table1_prefix + "_relations", gb, gen_degs);
}

int main_test1(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");

	db.execute_cmd("DELETE FROM E4b1_generators;");
	db.execute_cmd("DELETE FROM E4b1_relations;");
	db.execute_cmd("DELETE FROM E4b2_generators;");
	db.execute_cmd("DELETE FROM E4b2_relations;");
	db.execute_cmd("DELETE FROM E4b3_generators;");
	db.execute_cmd("DELETE FROM E4b3_relations;");
	db.execute_cmd("DELETE FROM E4b4_generators;");
	db.execute_cmd("DELETE FROM E4b4_relations;");
	db.execute_cmd("DELETE FROM E4b5_generators;");
	db.execute_cmd("DELETE FROM E4b5_relations;");
	db.execute_cmd("DELETE FROM E4b6_generators;");
	db.execute_cmd("DELETE FROM E4b6_relations;");

	auto start = std::chrono::system_clock::now();
	int t_max = 74;
	std::cout << "E4b1\n";
	generate_E4bk(db, "E4", "E4b1", { {{0, 1}} }, 58, t_max);
	std::cout << "E4b2\n";
	generate_E4bk(db, "E4b1", "E4b2", { {{54, 1}} }, 59, t_max);
	std::cout << "E4b3\n";
	generate_E4bk(db, "E4b2", "E4b3", { {{75, 1}} }, 60, t_max);
	std::cout << "E4b4\n";
	generate_E4bk(db, "E4b3", "E4b4", { {{113, 1}} }, 61, t_max);
	std::cout << "E4b5\n";
	generate_E4bk(db, "E4b4", "E4b5", { {{172, 1}} }, 62, t_max);
	std::cout << "E4b6\n";
	generate_E4bk(db, "E4b5", "E4b6", { {{253, 1}} }, 63, t_max);



	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Elapsed time: " << elapsed.count() << "s\n";

	return 0;
}

int main_generate_E4t(int argc, char** argv)
{
	return main_test1(argc, argv);

	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\ss.db)");

	auto start = std::chrono::system_clock::now();

	int t_max = 189;
	/*std::cout << "E4b1\n";
	generate_E4bk(conn, "E4", "E4b1", { {0, 1} }, 58, t_max);*/
	//std::cout << "E4b2\n";
	//generate_E4bk(conn, "E4b1", "E4b2", { {431, 1} }, 59, t_max);
	//std::cout << "E4b3\n";
	//generate_E4bk(conn, "E4b2", "E4b3", { {646, 1} }, 60, t_max);
	//std::cout << "E4b4\n";
	//generate_E4bk(conn, "E4b3", "E4b4", { {1064, 1} }, 61, t_max);
	std::cout << "E4b5\n";
	generate_E4bk(db, "E4b4", "E4b5", { {{1849, 1}} }, 62, t_max);
	/*std::cout << "E4b6\n";
	generate_E4bk(conn, "E4b5", "E4b6", { {252, 1} }, 63, t_max);*/
	//std::cout << "E4b7\n";
	//generate_E4bk(conn, "E4b6", "E4b7", 64);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Elapsed time: " << elapsed.count() << "s\n";

	return 0;
}