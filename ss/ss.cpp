#include "main.h"
#include "groebner.h"
#include "benchmark.h"
#include "linalg.h"
#include <iostream>

/* generate the basis of the algebra from the Groebner basis */
void generate_basis(const Database& db, const std::string& table_prefix, int t_max, bool drop_existing /*= false*/)
{
	/* load gen_degs, leadings, basis */
	if (drop_existing)
		db.execute_cmd("DELETE FROM " + table_prefix + "_basis;");
	std::vector<Deg> gen_degs = db.load_gen_degs(table_prefix + "_generators");
	Mon2d leadings = db.load_leading_terms(table_prefix + "_relations", t_max);
	std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", t_max);

	Timer timer;
	 
	/* starting t value */
	int t_min;
	if (!basis.empty())
		t_min = basis.rbegin()->first.t + 1;
	else {
		basis[Deg{ 0, 0, 0 }].push_back({}); /* if no monomial present insert the unit */
		t_min = 1;
		db.execute_cmd("INSERT INTO " + table_prefix + "_basis (mon_id, mon, s, t, v) VALUES (0, \"\", 0, 0, 0);");
	}

	/* Add new basis */
	for (int t = t_min; t <= t_max; t++) {
		std::map<Deg, Mon1d> basis_new;
		std::cout << "t=" << t << "          \r";
		for (int gen_id = (int)gen_degs.size() - 1; gen_id >= 0; --gen_id) {
			int t1 = t - gen_degs[gen_id].t;
			if (t1 >= 0) {
				auto p1 = basis.lower_bound(Deg{ 0, t1, 0 });
				auto p2 = basis.lower_bound(Deg{ 0, t1 + 1, 0 });
				for (auto p = p1; p != p2; ++p) {
					for (const auto& m : p->second) {
						if (m.empty() || gen_id <= m[0].gen) {
							Mon mon(mul(m, { {gen_id, 1} }));
							if ((size_t)gen_id >= leadings.size() || std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(),
								[&mon](const Mon& _m) { return divides(_m, mon); }))
								basis_new[p->first + gen_degs[gen_id]].push_back(std::move(mon));
						}
					}
				}
			}
		}
		/* Insert the new basis in degree t into the database */
		db.begin_transaction();
		db.save_basis(table_prefix + "_basis", basis_new);
		db.end_transaction();
		basis.merge(basis_new);
	}
}

/* generate the differetials on the monomials */
void generate_mon_diffs(const Database& db, const std::string& table_prefix, int r)
{
	/* load gen_degs, gen_diffs, gb and basis */

	std::vector<Deg> gen_degs = db.load_gen_degs(table_prefix + "_generators");
	Poly1d gen_diffs = db.load_gen_diffs(table_prefix + "_generators");
	Poly1d gb = db.load_gb(table_prefix + "_relations", -1);
	std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", -1);
	std::map<Deg, Poly1d> mon_diffs = db.load_mon_diffs(table_prefix + "_basis", basis, r, -1);

	/* compute diffs */

	Statement stmt(db, "UPDATE " + table_prefix + "_basis SET diff=?1 WHERE mon=?2;");

	db.begin_transaction();
	int prev_t = 0;
	for (auto& [d, basis_d] : basis) {
		if (mon_diffs.find(d) != mon_diffs.end())
			continue;
		if (d.t > prev_t) {
			prev_t = d.t;
			std::cout << "computing diffs, t=" << d.t << "          \r";
			db.end_transaction();
			db.begin_transaction();
		}
		for (size_t i = 0; i < basis_d.size(); ++i) {
			array diff_indices;
			if (basis_d[i].empty())
				mon_diffs[d].push_back({});
			else {
				int gen_id = basis_d[i].front().gen;
				Mon mon1 = div(basis_d[i], { {gen_id, 1} });
				Deg d1 = d - gen_degs[gen_id];
				size_t index_mon1 = std::lower_bound(basis.at(d1).begin(), basis.at(d1).end(), mon1) - basis.at(d1).begin();
				mon_diffs[d].push_back(grbn::Reduce(add(mul(gen_diffs[gen_id], mon1), mul({ {gen_id, 1} }, mon_diffs[d1][index_mon1])), gb));
				Deg d_diff = d + Deg{ 1, 0, -r };
				diff_indices = Poly_to_indices(mon_diffs[d][i], basis[d_diff]);
			}
			std::string str_diff(array_to_str(diff_indices));
			stmt.bind_str(1, str_diff);
			stmt.bind_str(2, Mon_to_str(basis_d[i]));
			stmt.step_and_reset();
		}
	}
	db.end_transaction();
}

/* generate the table of the spectral sequence */
void generate_ss(const Database& db, const std::string& table_basis, const std::string& table_ss, int r)
{
	std::map<Deg, array2d> mon_diffs_ind = db.load_mon_diffs_ind(table_basis, -1, true);
	std::map<Deg, Staircase> basis_ss;

	/* fill basis_ss */
	int prev_t = 0;
	for (auto& [deg, mon_diffs_d] : mon_diffs_ind) {
		if (deg.t > prev_t) {
			prev_t = deg.t;
			std::cout << "Compute basis_ss, t=" << deg.t << "          \r";
		}

		array all_indices = grbn::range((int)mon_diffs_d.size());
		array indices_with_diffs;

		array2d non_trivial_mon_diffs_d;
		for (int i = 0; i < (int)mon_diffs_d.size(); ++i)
			if (mon_diffs_d[i] != array{ -1 }) {
				non_trivial_mon_diffs_d.push_back(mon_diffs_d[i]);
				indices_with_diffs.push_back(i);
			}
		Staircase& basis_ss_d = basis_ss[deg];

		array lead_indices_of_boudaries;
		for (const array& base : basis_ss_d.basis_ind)
			if (base != array{ -1 })
				lead_indices_of_boudaries.push_back(base.front());
		std::sort(lead_indices_of_boudaries.begin(), lead_indices_of_boudaries.end());

		array indices_non_boundaries = lina::AddVectors(indices_with_diffs, lead_indices_of_boudaries);

		array2d map_diff;
		for (int i : indices_non_boundaries)
			map_diff.push_back(non_trivial_mon_diffs_d[i]);

		array2d image, kernel, g;
		lina::SetLinearMapV2(indices_non_boundaries, map_diff, image, kernel, g);

		/* fill Im(d_r) */
		Deg deg_diff = deg + Deg{ 1, 0, -r };
		for (size_t i = 0; i < image.size(); ++i) {
			basis_ss[deg_diff].basis_ind.push_back(std::move(image[i]));
			basis_ss[deg_diff].diffs_ind.push_back(std::move(g[i]));
			basis_ss[deg_diff].levels.push_back(r);
		}

		/* fill E_{r+2} */
		array lead_kernel;
		for (auto& cycle : kernel) {
			lead_kernel.push_back(cycle.front());
			basis_ss_d.basis_ind.push_back(cycle);
			basis_ss_d.diffs_ind.push_back({});
			basis_ss_d.levels.push_back(kLevelMax - r - 2);
		}
		std::sort(lead_kernel.begin(), lead_kernel.end());

		/* fill E_r */
		array rest_indices_with_diffs = lina::AddVectors(indices_non_boundaries, lead_kernel);
		for (int i : rest_indices_with_diffs) {
			basis_ss_d.basis_ind.push_back({ i });
			basis_ss_d.diffs_ind.push_back(std::move(mon_diffs_ind[deg][i]));
			basis_ss_d.levels.push_back(kLevelMax - r);
		}
		array indices_without_diffs = lina::AddVectors(all_indices, indices_with_diffs);
		for (int i : indices_without_diffs) {
			basis_ss_d.basis_ind.push_back({ i });
			basis_ss_d.diffs_ind.push_back({ -1 });
			basis_ss_d.levels.push_back(kLevelMax - r);
		}
	}

	basis_ss[Deg{ 0, 0, 0 }].diffs_ind = { {} };
	basis_ss[Deg{ 0, 0, 0 }].levels = { kLevelMax / 2 };

	/* insert into the database */
	db.begin_transaction();
	db.execute_cmd("DELETE FROM " + table_ss + ";");
	db.save_basis_ss(table_ss, basis_ss);
	db.end_transaction();
}

/* generate the homology of the E_r page */
void generate_next_page(const Database& db, const std::string& table_prefix, const std::string& table_H_prefix, int r)
{
	Poly1d gb = db.load_gb(table_prefix + "_relations", -1);
	std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", -1);
	std::map<Deg, BasisComplex> basis_ss = db.load_basis_ss(table_prefix + "_ss", r, -1);

	std::vector<Deg> gen_degs_H;
	array2d reprs_H;
	Poly1d gb_H;

	Mon2d leadings;
	std::map<Deg, Mon1d> basis_H;
	basis_H[Deg{ 0, 0, 0 }].push_back({});
	std::map<Deg, array2d> mon_reprs_H;
	mon_reprs_H[Deg{ 0, 0, 0 }].push_back({});
	std::map<Deg, Mon1d> basis_H_new;
	std::map<Deg, array2d> mon_reprs_H_new;

	int t_max = basis_ss.crbegin()->first.t;
	for (int t = 1; t <= t_max; ++t) {
		std::cout << "generate next page, t=" << t << "          \r";
		auto p1_ss = basis_ss.lower_bound(Deg{ 0, t, 0 });
		auto p2_ss = basis_ss.lower_bound(Deg{ 0, t + 1, 0 });
		for (auto p_ss = p1_ss; p_ss != p2_ss; ++p_ss) {
			const Deg& deg = p_ss->first;
			auto pNew = basis_H_new.find(deg);
			if (pNew != basis_H_new.end()) {
				/* Compute gb */
				array indices = grbn::range(int(pNew->second.size()));
				std::sort(indices.begin(), indices.end(), [&pNew](int a, int b) {return pNew->second[a] < pNew->second[b]; });
				Mon1d mons_new;
				array2d reprs_new;
				for (size_t j = 0; j < indices.size(); j++) {
					mons_new.push_back(std::move(pNew->second[indices[j]]));
					reprs_new.push_back(std::move(mon_reprs_H_new[deg][indices[j]]));
				}

				array2d image, kernel, g, quotient;
				lina::SetLinearMap(reprs_new, image, kernel, g);
				array lead_kernel;
				for (array& p_indices : kernel) {
					gb_H.push_back(indices_to_Poly(p_indices, mons_new));
					int index = gb_H.back().front().front().gen;
					if (size_t(index) >= leadings.size())
						leadings.resize(size_t(index) + 1);
					leadings[index].push_back(gb_H.back().front());
					lead_kernel.push_back(p_indices[0]);
				}
				std::sort(lead_kernel.begin(), lead_kernel.end());

				/* Find new generators and add to basis_H */
				quotient = lina::QuotientSpace(p_ss->second.cycles, image);
				lina::SimplifySpace(quotient);
				for (auto& x : quotient) {
					gen_degs_H.push_back(deg);
					reprs_H.push_back(x);
					basis_H[deg].push_back({ {int(gen_degs_H.size()) - 1, 1} });
					mon_reprs_H[deg].push_back(std::move(x));
				}

				/* Add to basis_H */
				array index_basis = lina::AddVectors(grbn::range(int(mons_new.size())), lead_kernel);
				for (int i : index_basis) {
					basis_H[deg].push_back(std::move(mons_new[i]));
					mon_reprs_H[deg].push_back(std::move(reprs_new[i]));
				}
			}
			else {
				/* Find new generators and add to basis_H */
				for (array& x : lina::SimplifySpace(p_ss->second.cycles)) {
					gen_degs_H.push_back(deg);
					reprs_H.push_back(x);
					basis_H[deg].push_back({ {int(gen_degs_H.size()) - 1, 1} });
					mon_reprs_H[deg].push_back(std::move(x));
				}
			}
		}
		for (auto& [deg, basis_H_new_d] : basis_H_new) {
			if (basis_ss.find(deg) == basis_ss.end()) {
				/* Add to gb_H */
				for (auto& mon : basis_H_new_d) {
					gb_H.push_back({ {std::move(mon)} });
					int index = gb_H.back().front().front().gen;
					if (size_t(index) >= leadings.size())
						leadings.resize(size_t(index) + 1);
					leadings[index].push_back(gb_H.back().front());
				}
			}
		}

		/* Compute new basis of degree t + 1 */
		if (t < t_max) {
			basis_H_new.clear();
			for (int gen_id = (int)gen_degs_H.size() - 1; gen_id >= 0; --gen_id) {
				int t1 = t + 1 - gen_degs_H[gen_id].t;
				if (t1 >= 0) {
					auto p1 = basis_H.lower_bound(Deg{ 0, t1, 0 });
					auto p2 = basis_H.lower_bound(Deg{ 0, t1 + 1, 0 });
					for (auto p = p1; p != p2; ++p) {
						const Deg& deg1 = p->first;
						const Deg& g_deg = gen_degs_H[gen_id];
						Deg deg = deg1 + g_deg;
						for (size_t j = 0; j < p->second.size(); ++j) {
							const auto& mon1 = p->second[j];
							const auto& repr1 = mon_reprs_H[deg1][j];
							if (mon1.empty() || gen_id <= mon1[0].gen) {
								Mon mon(mul(mon1, { {gen_id, 1} }));
								if (gen_id >= (int)leadings.size() || std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(),
									[mon](Mon _m) { return divides(_m, mon); })) {
									/* Compute the represeting cycle of the base monomial */
									auto poly1 = indices_to_Poly(repr1, basis[deg1]);
									auto gen_repr = indices_to_Poly(reprs_H[gen_id], basis[g_deg]);
									Poly repr = grbn::Reduce(mul(poly1, gen_repr), gb);
									array repr_indices = lina::Residue(basis_ss[deg].boundaries, Poly_to_indices(repr, basis[deg]));
									basis_H_new[deg].push_back(std::move(mon));
									mon_reprs_H_new[deg].push_back(std::move(repr_indices));
								}
							}
						}
					}
				}
			}
		}
	}

	db.begin_transaction();

	/* Save generators */
	Poly1d reprs_poly_H;
	for (size_t i = 0; i < (int)reprs_H.size(); ++i)
		reprs_poly_H.push_back(indices_to_Poly(reprs_H[i], basis[gen_degs_H[i]]));
	db.save_generators(table_H_prefix + "_generators", gen_degs_H, reprs_poly_H);

	/* Save gb */
	db.save_gb(table_H_prefix + "_relations", gb_H, gen_degs_H);

	/* Save basis */
	db.save_basis(table_H_prefix + "_basis", basis_H, mon_reprs_H);

	db.end_transaction();
}