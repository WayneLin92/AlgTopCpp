#include "main.h"

/* Add rels from gb in <= t to gb1 */
void AddRelFromGb(const Poly1d& gb, const array& gen_degs, Poly1d& gb1, RelHeap& heap1, int t, int t_max)
{
	/* Add relations from gb to gb1 */
	auto p1 = std::lower_bound(gb.begin(), gb.end(), t, [&gen_degs](const Poly& g, int t) {return get_deg(g, gen_degs) < t; });
	auto p2 = std::lower_bound(p1, gb.end(), t, [&gen_degs](const Poly& g, int t) {return get_deg(g, gen_degs) < t + 1; });
	for (auto pg = p1; pg != p2; ++pg)
		heap1.push(PolyWithT{ *pg, t });
	add_rels_from_heap(gb1, heap1, gen_degs, t, t_max);
}

/* compute ann * polys = 0 in degree t */
Poly2d ExtendAnn(const Poly1d& gb, const array& gen_degs, Poly1d& gb1, RelHeap& heap1, const Poly1d& polys, const array& deg_polys, int t, int t_max)
{
	/* Add relations from gb to gb1 */
	auto p1 = std::lower_bound(gb.begin(), gb.end(), t, [&gen_degs](const Poly& g, int t) {return get_deg(g, gen_degs) < t; });
	auto p2 = std::lower_bound(p1, gb.end(), t, [&gen_degs](const Poly& g, int t) {return get_deg(g, gen_degs) < t + 1; });
	for (auto pg = p1; pg != p2; ++pg)
		heap1.push(PolyWithT{ *pg, t });

	/* Add relations Xi=polys[i] to gb1 */
	int N = (int)polys.size();
	for (int i = N - 1; i >= 0; --i) {
		if (deg_polys[i] != t)
			break;
		Poly p = polys[i];
		p.push_back({ {-i - 1, 1} });
		heap1.push(PolyWithT{ std::move(p), t });
	}
	add_rels_from_heap(gb1, heap1, gen_degs, deg_polys, t, t_max);

	Poly2d result;
	if (polys.empty())
		return result;

	/* Extract linear relations from gb1 */
	for (auto pg = gb1.rbegin(); pg != gb1.rend(); ++pg) {
		if (get_deg(*pg, gen_degs) != t)
			break;
		if (pg->front()[0].gen < 0) {
			Poly1d ann;
			ann.resize(N);
			for (const Mon& m : *pg) {
				MonInd p = m.begin();
				for (; p != m.end() && p->gen < 0; ++p);
				Mon m1(m.begin(), p), m2(p, m.end());
				ann[size_t(-m1[0].gen) - 1] += reduce(mul(evaluate({ div(m1, { {m1[0].gen, 1} }) }, [&polys](int i) {return polys[size_t(-i) - 1]; }, gb), m2), gb);
			}
			result.push_back(std::move(ann));
		}
	}

	/* Add commutators to linear relations */
	for (int i = 0; i < N; ++i) {
		for (int j = i + 1; j < N; ++j) {
			if (deg_polys[i] + deg_polys[j] == t) {
				Poly1d result_i;
				result_i.resize(N);
				result_i[i] = polys[j];
				result_i[j] = polys[i];
				result.push_back(std::move(result_i));
			}
		}
	}

	return result;
}

/* Remove decomposables */
void Indecomposables(const Poly1d& gb, const array& gen_degs, Poly1d& gb1, RelHeap& heap1, Poly2d& vectors, const array& basis_degs, int t, int t_max)
{
	/* Add relations from gb to gb1 */
	auto p1 = std::lower_bound(gb.begin(), gb.end(), t, [&gen_degs](const Poly& g, int t) {return get_deg(g, gen_degs) < t; });
	auto p2 = std::lower_bound(p1, gb.end(), t, [&gen_degs](const Poly& g, int t) {return get_deg(g, gen_degs) < t + 1; });
	for (auto pg = p1; pg != p2; ++pg)
		heap1.push(PolyWithT{ *pg, t });

	/* Convert each vector v to a relation \\sum vi x_{-i-1} */
	Poly1d rels;
	for (const Poly1d& v : vectors) {
		Poly rel;
		for (int i = 0; i < basis_degs.size(); ++i)
			if (!v[i].empty())
				rel += v[i] * Mon{ {-i - 1, 1} };
		rels.push_back(std::move(rel));
	}
	add_rels_freemodule(gb1, heap1, FnGetDegV2{ gen_degs, basis_degs }, t, t_max);

	/* Add relations ordered by degree to gb1 */
	for (size_t i = 0; i < vectors.size(); ++i) {
		Poly rel = reduce(rels[i], gb1);
		if (!rel.empty())
			heap1.push(PolyWithT{ std::move(rel), t });
		else
			vectors[i].clear();
	}
	add_rels_freemodule(gb1, heap1, FnGetDegV2{ gen_degs, basis_degs }, t, t_max);

	/* Keep only the indecomposables in `vectors` */
	RemoveEmptyElements(vectors);
}

std::map<Deg, Mon1d> ExtendBasis(const std::vector<Deg>& gen_degs, const Mon2d& leadings, std::map<Deg, Mon1d>& basis, int t)
{
	std::map<Deg, Mon1d> basis_t;
	if (t == 0) {
		basis[Deg{ 0, 0, 0 }].push_back({});
		return basis_t;
	}
	for (int gen_id = (int)gen_degs.size() - 1; gen_id >= 0; --gen_id) {
		int t1 = t - gen_degs[gen_id].t;
		if (t1 >= 0) {
			auto p1_basis = basis.lower_bound(Deg{ 0, t1, 0 });
			auto p2_basis = basis.lower_bound(Deg{ 0, t1 + 1, 0 });
			for (auto p_basis = p1_basis; p_basis != p2_basis; ++p_basis) {
				for (auto p_m = p_basis->second.begin(); p_m != p_basis->second.end(); ++p_m) {
					if (p_m->empty() || gen_id <= p_m->front().gen) {
						Mon mon(mul(*p_m, { {gen_id, 1} }));
						if ((size_t)gen_id >= leadings.size() || std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(),
							[&mon](const Mon& _m) { return divides(_m, mon); }))
							basis_t[p_basis->first + gen_degs[gen_id]].push_back(std::move(mon));
					}
				}
			}
		}
	}
	return basis_t;
}

/* Compute relations in HB where B=A[X] in degrees (s, t, v) with fixed t and v2s=v+2*s */
std::vector<PolyWithT> FindRels(std::map<Deg, Mon1d>& basis_HB, const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X, const Poly1d& gb_A, const Poly1d& gen_reprs_B,
	int t, int v2s, const array& v_s)
{
	Mon1d basis_B_s; /* in deg s */
	Poly1d diffs_B_s; /* in deg s + 1 */
	Mon1d basis_B_sm1; /* in deg s - 1 */
	Poly1d diffs_B_sm1; /* in deg s */
	std::vector<PolyWithT> result;
	for (size_t i = 0; i < v_s.size(); ++i) {
		int s = v_s[i];
		Deg d = { s, t, v2s - 2 * s };
		if (i > 0 && v_s[i - 1] == s - 1) {
			std::swap(basis_B_sm1, basis_B_s); basis_B_s.clear();
			std::swap(diffs_B_sm1, diffs_B_s); diffs_B_s.clear();
		}
		else {
			get_basis_with_diff_B(basis_A, basis_X, basis_B_sm1, diffs_B_sm1, d - Deg{ 1, 0, -2 });
		}
		if (i < v_s.size() - 1 && v_s[i + 1] == s + 1)
			get_basis_with_diff_B(basis_A, basis_X, basis_B_s, diffs_B_s, d);
		else
			get_basis_B(basis_A, basis_X, basis_B_s, d);
		std::sort(basis_B_s.begin(), basis_B_s.end());

		array2d map_diff;
		for (Poly& diff : diffs_B_sm1)
			map_diff.push_back(Poly_to_indices(reduce(diff, gb_A), basis_B_s));
		array2d image_diff, kernel_diff, g_diff;
		set_linear_map(map_diff, image_diff, kernel_diff, g_diff);

		array2d map_repr;
		for (const Mon& mon : basis_HB[d]) {
			array repr = residue(image_diff, Poly_to_indices(get_repr({ mon }, gen_reprs_B, gb_A), basis_B_s));
			map_repr.push_back(std::move(repr));
		}
		array2d image_repr, kernel_repr, g_repr;
		set_linear_map(map_repr, image_repr, kernel_repr, g_repr);

		for (const array& rel_indices : kernel_repr)
			result.push_back(PolyWithT{ indices_to_Poly(rel_indices, basis_HB[d]), d.t });
		for (const array& rel_indices : kernel_repr)
			basis_HB[d][rel_indices[0]].clear();

		RemoveEmptyElements(basis_HB[d]);
	}
	return result;
}

void generate_HB(const Database& db, int t_max)
{
	/*# Load data */
	size_t index_x = (size_t)db.get_int("SELECT COUNT(*) FROM A_generators;"); /* the gen_id of x1 is index_x + 1 */
	size_t n = (size_t)db.get_int("SELECT COUNT(*) FROM B_generators;") - index_x;
	int t_min;
	try { t_min = db.get_int("SELECT MAX(t) FROM HA1_basis;") + 1; }
	catch (const char*) { t_min = 0; }
	std::vector<Deg> gen_degs_B = db.load_gen_degs("B_generators");
	Poly1d gen_diffs_B = db.load_gen_diffs("B_generators");
	Poly1d gb_A0 = db.load_gb("A_relations");
	std::map<Deg, DgaBasis1> basis_A0 = load_dga_basis(db, "A_basis", 2);
	std::vector<std::map<Deg, DgaBasis1>> basis_X;

	/* Data that needs to be updated */
	std::vector<std::vector<Deg>> gen_degs_HA;
	std::vector<array> gen_degs_t_HA; /* derived from gen_degs_HA */
	std::vector<Poly1d> gen_reprs_HA;
	std::vector<Poly1d> gb_HA;
	std::vector<Mon2d> leadings_HA; /* derived from gb_HA */
	std::vector<RelHeap> heap_HA; /* derived from gb_HA */
	std::vector<std::map<Deg, Mon1d>> basis_HA;
	std::vector<Poly1d> y;
	std::vector<array> t_y;

	std::vector<Poly1d> gb_HA_ann_c;
	std::vector<Poly1d> gb_HA_ind_y;
	std::vector<Poly1d> gb_HA_ann_y;
	std::vector<Poly1d> gb_HA_ind_a;
	std::vector<RelHeap> heap_HA_ann_c; /* derived from gb_HA_ann_c */
	std::vector<RelHeap> heap_HA_ind_y; /* derived from gb_HA_ind_y */
	std::vector<RelHeap> heap_HA_ann_y; /* derived from gb_HA_ann_y */
	std::vector<RelHeap> heap_HA_ind_a; /* derived from gb_HA_ind_a */
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

		/* Generate heap_HA */
		heap_HA.push_back(GenerateHeap(gb_HA.back(), gen_degs_t_HA.back(), {}, t_min, t_max));

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
			/* Load y and t_y */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_y  (y TEXT, t SMALLINT);"); }
			catch (const char*) {}
			y.push_back({}); t_y.push_back({});
			load_y(db, table_prefix + "_y", y.back(), t_y.back());

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

			/* Generate heaps */
			int t_x = gen_degs_B[index_x + i + 1].t; /* x = x_{i + 1} */
			heap_HA_ann_c.push_back(GenerateHeap(gb_HA_ann_c.back(), gen_degs_t_HA.back(), { gen_degs_B[index_x + i + 1].t }, t_min, t_max));
			heap_HA_ind_y.push_back(GenerateHeap(gb_HA_ind_y.back(), gen_degs_t_HA.back(), { gen_degs_B[index_x + i + 1].t }, t_min, t_max));
			heap_HA_ann_y.push_back(GenerateHeap(gb_HA_ann_y.back(), gen_degs_t_HA.back(), t_y[i], t_min - t_x, t_max - t_x));
			heap_HA_ind_a.push_back(GenerateHeap(gb_HA_ind_a.back(), gen_degs_t_HA.back(), t_y[i], t_min - t_x, t_max - t_x));
		}
	}

	/*# Compute */
	for (int t = t_min; t <= t_max; ++t) {
		db.begin_transaction();
		for (size_t i = 0; i < (size_t)n; ++i) {
			/*## compute ann(c) and ann(y) */
			/* Add x_{-1}=c to gb_HA_ann_c in order to compute ann_HA(c) */
			std::string table_prefix = i == 0 ? "HA" : "HA" + std::to_string(i);
			Poly c; /* dx = c */
			Deg deg_x = gen_degs_B[index_x + i + 1] + Deg{ 1, 0, -2 };
			int t_x = deg_x.t; /* x = x_{i + 1} */
			if (t == t_x) {
				c = proj(gen_diffs_B[index_x + i + 1], gen_degs_B, gen_diffs_B, gb_A0, basis_A0, basis_X[i], gen_reprs_HA[i], basis_HA[i]);
				heap_HA[i + 1].push(PolyWithT{ c, t_x }); /* Add relation dx=0 */
			}

			Poly2d y_new; /* Annilators of c in gb_HA[i] in degree t - t_x */
			std::map<Deg, Mon1d> basis_t;
			if (t < t_x) {
				AddRelFromGb(gb_HA[i], gen_degs_t_HA[i], gb_HA_ann_c[i], heap_HA_ann_c[i], t, t_max);
				AddRelFromGb(gb_HA[i], gen_degs_t_HA[i], gb_HA_ann_y[i], heap_HA_ann_y[i], t, t_max);
				AddRelFromGb(gb_HA[i], gen_degs_t_HA[i], gb_HA[i + 1], heap_HA[i + 1], t, t_max);
				
				auto p1_basis_HA = basis_HA[i].lower_bound(Deg{ 0, t, 0 });
				auto p2_basis_HA = basis_HA[i].lower_bound(Deg{ 0, t + 1, 0 });
				for (auto p_basis_HA = p1_basis_HA; p_basis_HA != p2_basis_HA; ++p_basis_HA)
					basis_t[p_basis_HA->first] = p_basis_HA->second;
			}
			else {
				y_new = ExtendAnn(gb_HA[i], gen_degs_t_HA[i], gb_HA_ann_c[i], heap_HA_ann_c[i], { c }, { t_x }, t, t_max);
				Indecomposables(gb_HA[i], gen_degs_t_HA[i], gb_HA_ind_y[i], heap_HA_ind_y[i], y_new, { t_x }, t, t_max);
				for (auto& yi : y_new) {
					y[i].push_back(yi[0]);
					t_y[i].push_back(t - t_x);
				}
				save_y(db, table_prefix + "_y", y_new, t);
				Poly2d a = ExtendAnn(gb_HA[i], gen_degs_t_HA[i], gb_HA_ann_c[i], heap_HA_ann_c[i], y[i], t_y[i], t - t_x, t_max - t_x);
				Indecomposables(gb_HA[i], gen_degs_t_HA[i], gb_HA_ind_a[i], heap_HA_ind_a[i], a, t_y[i], t - t_x, t_max - t_x);

				/* Add new generators to HA[i+1] */
				Poly1d cy_inv;
				for (const Poly1d& yk : y_new) {
					Poly cyk = c * yk[0];
					Poly cyk_repr = get_repr(cyk, gen_reprs_HA[i], gb_A0);
					cy_inv.push_back(d_inv(cyk_repr, gen_degs_B, gen_diffs_B, gb_A0, basis_A0, basis_X[i]));
				}
				int gen_degs_HA_t_start = (int)gen_degs_HA[i + 1].size();
				if (t == t_x * 2) { /* Add [x^2] */
					gen_degs_HA[i + 1].push_back(deg_x * 2);
					gen_degs_t_HA[i + 1].push_back(gen_degs_HA[i + 1].back().t);
					gen_reprs_HA[i + 1].push_back({ {{int(index_x + i + 1), 2}} });
				}
				for (size_t j = 0; j < y_new.size(); ++j) { /* Add [xy+d^{-1}(cy)] */
					gen_degs_HA[i + 1].push_back(deg_x + get_deg(y_new[i][j], gen_degs_HA[i]));
					gen_degs_t_HA[i + 1].push_back(gen_degs_HA[i + 1].back().t);
					gen_reprs_HA[i + 1].push_back(Mon{ {int(index_x + i + 1), 2} } *get_repr(y_new[i][j], gen_reprs_HA[i], gb_A0) + cy_inv[j]);
				}
				/* Add new relations to HA[i + 1] */
				/* Degrees of relations */
				std::map<int, array> rel_degs;
				for (size_t j = gen_degs_HA[i + 1].size() - y.size(); j < gen_degs_HA[i + 1].size(); ++j) // Needs optimization
					for (size_t k = j; k < gen_degs_HA[i + 1].size(); ++k) {
						Deg deg_gjgk = gen_degs_HA[i + 1][j] + gen_degs_HA[i + 1][k];
						if (deg_gjgk.t == t)
							rel_degs[deg_gjgk.v + 2 * deg_gjgk.s].push_back(deg_gjgk.s);
					}
				for (const Poly1d& ak : a)
					for (size_t j = 0; j < y.size(); ++j)
						if (!ak[j].empty()) {
							Deg deg_akgj = get_deg(ak[j], gen_degs_HA[i + 1]) + gen_degs_HA[i + 1][gen_degs_HA[i + 1].size() - y.size() + j];
							rel_degs[deg_akgj.v + 2 * deg_akgj.s].push_back(deg_akgj.s);
							break;
						}
				for (auto& [v1, v_s] : rel_degs) {
					std::sort(v_s.begin(), v_s.end());
					v_s.erase(std::unique(v_s.begin(), v_s.end()), v_s.end());
				}

				/* Compute relations */
				leadings_HA[i + 1].clear();
				leadings_HA[i + 1].resize(gen_degs_HA[i + 1].size());
				for (const Poly& g : gb_HA[i + 1])
					leadings_HA[i + 1][g[0][0].gen].push_back(g[0]);
				basis_t = ExtendBasis(gen_degs_HA[i + 1], leadings_HA[i + 1], basis_HA[i + 1], t);
				for (auto& [v2s, v_s] : rel_degs) {
					std::vector<PolyWithT> rels = FindRels(basis_t, basis_A0, basis_X[i + 1], gb_A0, gen_reprs_HA[i + 1], t, v2s, v_s);
					for (PolyWithT& rel : rels)
						heap_HA[i + 1].push(rel);
				}
				add_rels_from_heap(gb_HA[i + 1], heap_HA[i + 1], gen_degs_t_HA[i + 1], t, t_max);
			}

			/* Save data */
			db.save_generators(table_prefix + "_generators", gen_degs_HA[i + 1], gen_reprs_HA[i + 1]);
			db.save_gb(table_prefix + "_gb", gb_HA[i + 1], gen_degs_HA[i + 1]);
			db.save_basis(table_prefix + "_basis", basis_t);
			basis_HA[i + 1].merge(basis_t);
			SaveGb(db, table_prefix + "_ann_c_gb", gb_HA_ann_c[i], gen_degs_t_HA[i], { t_x }, t);
			SaveGb(db, table_prefix + "_ind_y_gb", gb_HA_ann_c[i], gen_degs_t_HA[i], { t_x }, t);
			SaveGb(db, table_prefix + "_ann_y_gb", gb_HA_ann_c[i], gen_degs_t_HA[i], t_y[i], t - t_x);
			SaveGb(db, table_prefix + "_ind_a_gb", gb_HA_ann_c[i], gen_degs_t_HA[i], t_y[i], t - t_x);

			std::vector<Poly1d> gb_HA_ann_c;
			std::vector<Poly1d> gb_HA_ind_y;
			std::vector<Poly1d> gb_HA_ann_y;
			std::vector<Poly1d> gb_HA_ind_a;
		}
		db.end_transaction();
	}


}

int main_test1(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\HB.db)");

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