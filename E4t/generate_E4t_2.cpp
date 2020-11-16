#include "../algtop/main.h"
#include <fstream>
#include <future>

/********** STRUCTS AND CLASSES **********/

struct DgaBasis1
{
	Mon1d basis;
	Poly1d diffs;
};

/********** FUNCTIONS **********/

inline Poly get_repr(const Mon& mon, const Poly1d& gen_reprs, const Poly1d& gb)
{
	return evaluate({ mon }, [&gen_reprs](int i) {return gen_reprs[i]; }, gb);
}

void ExtendAnn(const Poly1d& gb, const Poly1d& polys, const array& gen_degs, int t, int deg_max) {};//

/* return a basis for polynomails of xi for t<=t_max */
std::map<Deg, DgaBasis1> get_basis_X(const std::vector<Deg>& gen_degs, const Poly1d& gen_diffs, int start_x, int num_x, int t_max)
{
	std::map<Deg, DgaBasis1> result;
	result[Deg{ 0, 0, 0 }].basis.push_back({});
	result[Deg{ 0, 0, 0 }].diffs.push_back({});

	for (int t = 1; t <= t_max; ++t) {
		std::map<Deg, DgaBasis1> basis_new;
		std::cout << "Computing basis_x, t=" << t << "          \r";
		for (int gen_id = start_x + num_x - 1; gen_id >= start_x; --gen_id) { /* 58 is the gen_id of b1 */
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

void save_basis_X(Database& db, int start_x, int num_x)
{
	int t_max = 200;
	std::vector<Deg> gen_degs_B = db.load_gen_degs("B_generators");
	Poly1d diffs_B = db.load_gen_diffs("B_generators");
	std::map<Deg, DgaBasis1> basis_X = get_basis_X(gen_degs_B, diffs_B, start_x, num_x, t_max);

	std::string table_name = "X" + std::to_string(num_x) + "_basis";
	try { db.execute_cmd("DROP TABLE " + table_name + ';'); }
	catch (const char*) {}
	db.execute_cmd("CREATE TABLE " + table_name + " (mon_id INTEGER PRIMARY KEY, mon TEXT NOT NULL UNIQUE, diff TEXT, s SMALLINT, t SMALLINT, v SMALLINT);");
	Statement stmt;
	stmt.init(db, "INSERT INTO " + table_name + " (s, t, v, mon, diff) VALUES (?1, ?2, ?3, ?4, ?5);");

	db.begin_transaction();
	for (auto& [deg, basis_d] : basis_X) {
		for (size_t i = 0; i < basis_d.basis.size(); ++i) {
			stmt.bind_int(1, deg.s);
			stmt.bind_int(2, deg.t);
			stmt.bind_int(3, deg.v);
			stmt.bind_str(4, Mon_to_str(basis_d.basis[i]));
			stmt.bind_str(5, Poly_to_str(basis_X[deg].diffs[i]));
			stmt.step_and_reset();
		}
	}
	db.end_transaction();
	std::cout << table_name << " is created, number of degrees=" << basis_X.size() << '\n';
}


std::map<Deg, DgaBasis1> load_dga_basis1(const Database& db, const std::string& table_name, int r)
{
	std::map<Deg, DgaBasis1> basis;
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
	std::cout << "dga_basis loaded from " << table_name << ", size=" << basis.size() << '\n';
	return basis;
}

std::map<Deg, DgaBasis1> load_basis_X(const Database& db, const std::string& table_name, int t_max, int r)
{
	std::map<Deg, DgaBasis1> basis;
	Statement stmt;
	stmt.init(db, "SELECT mon, diff, s, t, v FROM " + table_name + " WHERE t<=" + std::to_string(t_max) + " ORDER BY mon_id;");
	int prev_t = 0;
	while (stmt.step() == SQLITE_ROW) {
		Deg d = { stmt.column_int(2), stmt.column_int(3), stmt.column_int(4) };
		basis[d].basis.push_back(str_to_Mon(stmt.column_str(0)));
		basis[d].diffs.push_back(str_to_Poly(stmt.column_str(1)));
	}
	std::cout << "basis loaded from " << table_name << ", size=" << basis.size() << '\n';
	return basis;
}

void get_basis_B(const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X, Mon1d* p_basis, Poly1d* p_diffs, const Deg& deg)
{
	for (auto p2 = basis_X.rbegin(); p2 != basis_X.rend(); ++p2) {
		const Deg& d2 = p2->first;
		if (d2.s <= deg.s && d2.t <= deg.t && d2.v <= deg.v) {
			const Deg d1 = deg - d2;
			auto p1 = basis_A.find(d1);
			if (p1 != basis_A.end()) {
				if (p_diffs == nullptr) {
					for (size_t i = 0; i < p1->second.basis.size(); ++i)
						for (size_t j = 0; j < p2->second.basis.size(); ++j)
							p_basis->push_back(mul(p1->second.basis[i], p2->second.basis[j]));
				}
				else {
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
}

/* Assume poly is a boundary. Return the chain with it as boundary */
Poly d_inv(const Poly& poly, const std::vector<Deg>& gen_degs, const Poly1d& diffs, const Poly1d& gb, const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X)
{
	if (poly.empty())
		return {};
	Deg deg_poly = get_deg(poly[0], gen_degs);
	Deg deg_result = deg_poly - Deg{ 1, 0, -2 };
	Mon1d basis_in_poly;
	Mon1d basis_in_result;
	get_basis_B(basis_A, basis_X, &basis_in_poly, nullptr, deg_poly);
	get_basis_B(basis_A, basis_X, &basis_in_result, nullptr, deg_result);
	std::sort(basis_in_poly.begin(), basis_in_poly.end());
	std::sort(basis_in_result.begin(), basis_in_result.end());
	array2d map_diff;
	for (const Mon& mon : basis_in_result)
		map_diff.push_back(Poly_to_indices(reduce(get_diff(mon, diffs), gb), basis_in_poly));
	array2d image, kernel, g;
	set_linear_map(map_diff, image, kernel, g);
	return indices_to_Poly(get_image(image, g, Poly_to_indices(poly, basis_in_poly)), basis_in_result);
}

/* Assume poly is a cycle. Return the homology class */
Poly proj(const Poly& poly, const std::vector<Deg>& gen_degs, const Poly1d& gen_diffs, const Poly1d& gb, const std::map<Deg, DgaBasis1>& basis_A,
	const std::map<Deg, DgaBasis1>& basis_X, const Poly1d& gen_reprs, std::map<Deg, Mon1d>& basis_H)
{
	if (poly.empty())
		return {};
	Deg deg_poly = get_deg(poly[0], gen_degs);
	Deg deg_result = deg_poly - Deg{ 1, 0, -2 };
	Mon1d basis_in_poly;
	Mon1d basis_in_result;
	get_basis_B(basis_A, basis_X, &basis_in_poly, nullptr, deg_poly);
	get_basis_B(basis_A, basis_X, &basis_in_result, nullptr, deg_result);
	std::sort(basis_in_poly.begin(), basis_in_poly.end());
	std::sort(basis_in_result.begin(), basis_in_result.end());

	array2d map_diff;
	for (const Mon& mon : basis_in_result)
		map_diff.push_back(Poly_to_indices(reduce(get_diff(mon, gen_diffs), gb), basis_in_poly));
	array2d image, kernel, g;
	set_linear_map(map_diff, image, kernel, g);

	array2d map_repr;
	for (const Mon& mon : basis_H[deg_poly])
		map_repr.push_back(residue(image, Poly_to_indices(get_repr(mon, gen_reprs, gb), basis_in_poly)));
	array2d image1, kernel1, g1;
	set_linear_map(map_repr, image1, kernel1, g1);

	return indices_to_Poly(get_image(image1, g1, Poly_to_indices(poly, basis_in_poly)), basis_H[deg_poly]);
}

//std::vector<PolyWithT> find_relations(Deg d, Mon1d& basis_d, const std::map<Deg, DgaBasis1>& basis_E2, const std::map<Deg, DgaBasis1>& basis_bi, const Poly1d& gb_E2t, const Poly1d& reprs)
//{
//	Mon1d basis_d_E2t;
//	get_basis_E2t(basis_E2, basis_bi, &basis_d_E2t, nullptr, d);
//	Mon1d basis_d1_E2t;
//	Poly1d diffs_d_E2t;
//	get_basis_E2t(basis_E2, basis_bi, &basis_d1_E2t, &diffs_d_E2t, d - Deg{ 1, 0, -2 });
//	array indices = range((int)basis_d1_E2t.size());
//	std::sort(indices.begin(), indices.end(), [&basis_d1_E2t](int i1, int i2) {return basis_d1_E2t[i1] < basis_d1_E2t[i2]; });
//	std::sort(basis_d_E2t.begin(), basis_d_E2t.end());
//	std::sort(basis_d1_E2t.begin(), basis_d1_E2t.end());
//
//	array2d map_diff;
//	for (int i : indices)
//		map_diff.push_back(Poly_to_indices(reduce(diffs_d_E2t[indices[i]], gb_E2t), basis_d_E2t));
//	array2d image_diff, kernel_diff, g_diff;
//	set_linear_map(map_diff, image_diff, kernel_diff, g_diff);
//
//	array2d map_repr;
//	for (const Mon& mon : basis_d) {
//		array repr = residue(image_diff, Poly_to_indices(evaluate({ mon }, [&reprs](int i) {return reprs[i]; }, gb_E2t), basis_d_E2t));
//		map_repr.push_back(std::move(repr));
//	}
//	array2d image_repr, kernel_repr, g_repr;
//	set_linear_map(map_repr, image_repr, kernel_repr, g_repr);
//
//	std::vector<PolyWithT> result;
//	for (const array& rel_indices : kernel_repr)
//		result.push_back(PolyWithT{ indices_to_Poly(rel_indices, basis_d), d.t });
//
//	for (const array& rel_indices : kernel_repr)
//		basis_d[rel_indices[0]].clear();
//	RemoveEmptyElements(basis_d);
//
//	return result;
//}

void generate_HB(const Database& db, int t_max)
{
	/*# Load data */
	int index_x = db.get_int("SELECT COUNT(*) FROM A_generators;"); /* the gen_id of x1 is index_x + 1 */
	int n = db.get_int("SELECT COUNT(*) FROM B_generators;") - index_x;
	int t_min;
	try { t_min = db.get_int("SELECT MAX(t) FROM HA1_basis;") + 1; }
	catch (const char*) { t_min = 0; }
	std::vector<Deg> gen_degs_B = db.load_gen_degs("B_generators");
	Poly1d gen_diffs_B = db.load_gen_diffs("B_generators");
	Poly1d gb_A = db.load_gb("A_relations");
	std::map<Deg, DgaBasis1> basis_A0 = load_dga_basis1(db, "A_basis", 2);
	std::vector<std::map<Deg, DgaBasis1>> basis_X;

	std::vector<std::vector<Deg>> gen_degs_HA;
	std::vector<array> gen_degs_t_HA;
	std::vector<Poly1d> gen_reprs_HA;
	std::vector<Poly1d> gb_HA;
	std::vector<Mon2d> leadings_HA;
	std::vector<std::map<Deg, Mon1d>> basis_HA;
	std::vector<array> t_y;

	std::vector<Poly1d> gb_HA_ann_c;
	std::vector<Poly1d> gb_HA_ann_y;
	std::vector<Poly1d> gb_HA_ind_y;
	std::vector<Poly1d> gb_HA_ind_a;
	std::vector<RelHeap> heap_HA_ann_c;
	std::vector<RelHeap> heap_HA_ann_y;
	std::vector<RelHeap> heap_HA_ind_y;
	std::vector<RelHeap> heap_HA_ind_a;
	for (int i = 0; i <= n; ++i) {
		std::string table_prefix = i == 0 ? "HA" : "HA" + std::to_string(i);

		/* Load gen_degs_HA */
		try { db.execute_cmd("CREATE TABLE " + table_prefix + "_generators (gen_id INTEGER PRIMARY KEY, gen_name TEXT UNIQUE, gen_diff TEXT, repr TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
		catch (const char*) {}
		gen_degs_HA.push_back(db.load_gen_degs(table_prefix + "_generators"));

		/* Initialize gen_degs_t_HA */
		gen_degs_t_HA.push_back({});
		for (auto p = gen_degs_HA.back().begin(); p < gen_degs_HA.back().end(); ++p)
			gen_degs_t_HA.back().push_back(p->t);

		/* Load gen_reprs_HA */
		gen_reprs_HA.push_back(db.load_gen_reprs(table_prefix + "_generators"));

		/* Load gb_HA */
		try { db.execute_cmd("CREATE TABLE " + table_prefix + "_relations (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
		catch (const char*) {}
		gb_HA.push_back(db.load_gb(table_prefix + "_relations"));

		/* Initialize leadings_HA */
		leadings_HA.push_back({});
		leadings_HA.back().resize(gen_degs_HA.back().size());
		for (const Poly& g : gb_HA.back())
			leadings_HA.back()[g[0][0].gen].push_back(g[0]);

		/* Load cache data for Groebner basis */
		if (i == 0) {
			basis_X.push_back({});
			basis_HA.push_back({});
		}
		else {
			/* Load basis_X */
			basis_X.push_back(load_basis_X(db, "X" + std::to_string(i) + "_basis", t_max, 2));

			/* Load basis_HA */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_basis  (mon_id INTEGER PRIMARY KEY, mon TEXT NOT NULL UNIQUE, diff TEXT, repr TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
			catch (const char*) {}
			basis_HA.push_back(db.load_basis(table_prefix + "_basis"));
		}
		if (i < n){
			/* Load t_y */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_t_y  (t SMALLINT);"); }
			catch (const char*) {}
			t_y.push_back(db.get_ints(table_prefix + "_t_y", "t"));

			/* Load gb_HA_ann_c */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_ann_c_gb (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
			catch (const char*) {}
			gb_HA_ann_c.push_back(db.load_gb(table_prefix + "_ann_c_gb"));

			/* Load gb_HA_ann_y */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_ann_y_gb (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
			catch (const char*) {}
			gb_HA_ann_y.push_back(db.load_gb(table_prefix + "_ann_y_gb"));

			/* Load gb_HA_ind_y */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_ind_y_gb (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
			catch (const char*) {}
			gb_HA_ind_y.push_back(db.load_gb(table_prefix + "_ind_y_gb"));

			/* Load gb_HA_ind_a */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_ind_a_gb (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
			catch (const char*) {}
			gb_HA_ind_a.push_back(db.load_gb(table_prefix + "_ind_a_gb"));

			/* Load heaps */
			heap_HA_ann_c.push_back(GenerateHeap(gb_HA_ann_c.back(), gen_degs_t_HA.back(), { gen_degs_B[index_x + i + 1].t }, t_min, t_max));
			heap_HA_ann_y.push_back(GenerateHeap(gb_HA_ann_y.back(), gen_degs_t_HA.back(), t_y[i], t_min, t_max));
			heap_HA_ind_y.push_back(GenerateHeap(gb_HA_ind_y.back(), gen_degs_t_HA.back(), { gen_degs_B[index_x + i + 1].t }, t_min, t_max));
			heap_HA_ind_a.push_back(GenerateHeap(gb_HA_ind_a.back(), gen_degs_t_HA.back(), t_y[i], t_min, t_max));
		}
	}

	/*# Compute Ann(c) and Ann(y) */
	for (size_t i = 0; i < (size_t)n; ++i) {
		int t_x = gen_degs_B[index_x + i + 1].t; /* x = x_{i + 1} */
		Poly c; /* dx = c */
		for (int t = t_min; t <= t_max; ++t) {
			/*## compute y */
			if (t == t_x) {
				/* Add x_{-1}=c to gb_HA_ann_c in order to compute ann_HA(c) */
				c = proj(gen_diffs_B[index_x + i], gen_degs_B, gen_diffs_B, gb_A, basis_A0, basis_X[i], gen_reprs_HA[i], basis_HA[i]);
				Poly p = add(c, Poly{ {{-1, 1}} });
				heap_HA_ann_c[i].push(PolyWithT{ p, t });
			}
			add_rels_from_heap(gb_HA_ann_c[i], heap_HA_ann_c[i], gen_degs_t_HA[i], array{ t_x }, t, t_max);

			Poly1d y; /* Annilators of c in gb_HA[i] in degree t - t_x */
			if (t > t_x) {
				for (auto pg = gb_HA_ann_c[i].rbegin(); pg != gb_HA_ann_c[i].rend(); ++pg) {
					if (get_deg(*pg, gen_degs_HA[i]).t != t)
						break;
					if (pg->front().front().gen < 0) {
						Poly ann;
						for (const Mon& m : *pg) {
							MonInd p = m.begin();
							for (; p != m.end() && p->gen < 0; ++p);
							Mon m1(m.begin(), p), m2(p, m.end());
							ann += reduce(evaluate({ div(m1, { {m1[0].gen, 1} }) }, [&c](int) {return c; }, gb_HA[i]) * m2, gb_HA[i]);
						}
						y.push_back(std::move(ann));
					}
				}
			}

			/* Remove yi which are decomposables */
			add_rels_from_heap(gb_HA_ind_y[i], heap_HA_ind_y[i], gen_degs_t_HA[i], array{ t_x }, t, t_max);
			for (auto&& yk : y) {
				Poly rel = reduce(yk, gb_HA_ind_y[i]);
				if (!rel.empty())
					heap_HA_ind_y[i].push(PolyWithT{ std::move(rel), t });
				else
					yk.clear();
			}
			RemoveEmptyElements(y);


		}
	}

	///* compute new generators */
	//Poly2d b_ = ann_seq(gb, { a }, gen_degs_t, t_max);
	//std::cout << "b_=" << b_ << '\n';
	//std::cout << "b_.size()=" << b_.size() << '\n';

	//Poly1d b;
	//for (Poly1d& v : b_)
	//	b.push_back(std::move(v[0]));
	//Poly2d c_ = ann_seq(gb, b, gen_degs_t, t_max - gen_degs_E2t[bk_gen_id].t);
	//std::cout << "c_.size()=" << c_.size() << '\n';

	//std::vector<Deg> gen_degs_E2bk_1(gen_degs_E2t);
	//gen_degs_E2bk_1.pop_back();
	//Poly1d ab_inv;
	//for (const Poly& bi : b) {
	//	Poly abi = mul(a, bi);
	//	Poly abi_repr = evaluate(abi, [&reprs](int i) {return reprs[i]; }, gb_E2t);
	//	ab_inv.push_back(d_inv(abi_repr, gb_E2t, leadings_E2t, gen_degs_E2bk_1, diffs_E2t));//
	//}

	//if (gen_degs_E2t.back().t * 2 <= t_max) { /* Add [x^2] */
	//	gen_degs.push_back(gen_degs_E2t.back() * 2);
	//	gen_degs_t.push_back(gen_degs.back().t);
	//	reprs.push_back({ {{bk_gen_id, 2}} });
	//}
	//for (size_t i = 0; i < b.size(); ++i) { /* Add [xb+d^{-1}(ab)] */
	//	gen_degs.push_back(gen_degs_E2t[bk_gen_id] + get_deg(b[i], gen_degs));
	//	gen_degs_t.push_back(gen_degs.back().t);
	//	reprs.push_back(add(mul(evaluate(b[i], [&reprs](int i) {return reprs[i]; }, gb_E2t), { {bk_gen_id, 1} }), ab_inv[i]));
	//}

	///* compute the relations */
	//add_rels_from_heap(gb, { a }, gen_degs_t, t_max); /* Add relations dx=0 */
	//leadings.clear();
	//leadings.resize(gen_degs.size());
	//for (const Poly& g : gb)
	//	leadings[g[0][0].gen].push_back(g[0]);

	//std::vector<Deg> rel_degs;
	//for (size_t i = gen_degs.size() - b.size(); i < gen_degs.size(); ++i)
	//	for (size_t j = i; j < gen_degs.size(); ++j) {
	//		Deg deg_gigj = gen_degs[i] + gen_degs[j];
	//		if (deg_gigj.t <= t_max)
	//			rel_degs.push_back(deg_gigj);
	//	}
	//for (const Poly1d& ci : c_)
	//	for (size_t j = 0; j < b.size(); ++j)
	//		if (!ci[j].empty()) {
	//			rel_degs.push_back(gen_degs[gen_degs.size() - b.size() + j] + get_deg(ci[j], gen_degs));
	//			break;
	//		}
	//std::sort(rel_degs.begin(), rel_degs.end());
	//rel_degs.erase(std::unique(rel_degs.begin(), rel_degs.end()), rel_degs.end());

	//std::ofstream myfile;
	//myfile.open("rel_gens.txt");
	//for (const auto& d : rel_degs)
	//	myfile << d << '\n';
	//myfile.close();

	//std::map<Deg, Mon1d> basis;
	//size_t i = 0;
	//std::vector<PolyWithT> heap;
	//for (int t = 1; t <= t_max; ++t) {
	//	get_basis(leadings, gen_degs, basis, t);

	//	std::vector<std::future<std::vector<PolyWithT>>> futures;
	//	for (; i < rel_degs.size() && rel_degs[i].t == t; ++i) {
	//		const Deg d = rel_degs[i];
	//		Mon1d& basis_d = basis[d];
	//		futures.push_back(std::async(std::launch::async, find_relations, d, std::ref(basis_d), std::ref(basis_A), std::ref(basis_X), std::ref(gb_E2t), std::ref(reprs)));
	//	}
	//	for (size_t j = 0; j < futures.size(); ++j) {
	//		futures[j].wait();
	//		std::cout << "t=" << t << " completed thread=" << j + 1 << '/' << futures.size() << "          \r";
	//	}
	//	for (auto& f : futures)
	//		for (auto& rel : f.get()) {
	//			heap.push_back(std::move(rel));
	//			std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
	//		}
	//	
	//	add_rels_from_heap(gb, heap, gen_degs_t, std::min(t + 1, t_max), t_max);
	//	leadings.clear();
	//	leadings.resize(gen_degs.size());
	//	for (const Poly& g : gb)
	//		leadings[g[0][0].gen].push_back(g[0]);
	//}

	///* Save generators and relations */
	//db.save_generators(table1_prefix + "_generators", gen_degs, reprs);
	//db.save_gb(table1_prefix + "_relations", gb, gen_degs);
}

int main_test1(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\HB.db)");

	/*for (int n = 1; n <= 7; n++) // Create Xi_basis
		save_basis_X(db, 58, n);
	return 0;*/

	return 0;
}

int main_generate_E4t(int argc, char** argv)
{
	//return main_test1(argc, argv);

	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\HB.db)");
	Timer timer;

	generate_HB(db, 74);

	return 0;
}