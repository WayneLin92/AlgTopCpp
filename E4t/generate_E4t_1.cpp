#include "main.h"
#include <chrono>
#include <future>
#include <set>

/********** STRUCTS AND CLASSES **********/

struct DgaBasis1
{
	Mon1d basis;
	Poly1d diffs;
};

/********** FUNCTIONS **********/

std::map<Deg, DgaBasis1> load_dga_basis(const Database& db, const std::string& table_name, int r)
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
	std::cout << "basis loaded from " << table_name << ", size=" << basis.size() << '\n';
	return basis;
}

std::map<Deg, DgaBasis1> load_basis_bi(const Database& db, const std::string& table_name, int t_max, int r)
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


void get_basis_E2t_1(const std::map<Deg, DgaBasis1>& basis_E2, const std::map<Deg, DgaBasis1>& basis_bi, Mon1d* p_basis, Poly1d* p_diffs, const Deg& deg)
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

void generate_E4t(const Database& db, int t_max)
{
	const std::string table_prefix = "E2t";
	const std::string table_H_prefix = "E4t";
	const int r = 2;
	Poly1d gb = db.load_gb(table_prefix + "_relations");
	std::map<Deg, DgaBasis1> basis_E2 = load_dga_basis(db, "E2_basis", r);
	std::map<Deg, DgaBasis1> basis_bi = load_basis_bi(db, "E2t_bi_basis", t_max, r);


	std::vector<Deg> gen_degs_H;
	Poly1d reprs_H;
	Poly1d gb_H;

	Mon2d leadings;
	std::map<Deg, Mon1d> basis_H;
	basis_H[Deg{ 0, 0, 0 }].push_back({});
	std::map<Deg, Poly1d> mon_reprs_H;
	mon_reprs_H[Deg{ 0, 0, 0 }].push_back({});
	std::map<Deg, Mon1d> basis_H_new;
	std::map<Deg, Poly1d> mon_reprs_H_new;

	for (int t = 1; t <= t_max; ++t) {
		std::cout << "generate next page, t=" << t << "          \n";

		std::set<Deg> degs;
		auto pE2End = basis_E2.lower_bound(Deg{ 0, t + 1, 0 });
		for (auto pE2 = basis_E2.begin(); pE2 != pE2End; ++pE2) {
			auto pBi = basis_bi.lower_bound(Deg{ 0, t - pE2->first.t, 0 });
			auto pBiEnd = basis_bi.lower_bound(Deg{ 0, t - pE2->first.t + 1, 0 });
			for (; pBi != pBiEnd; ++pBi)
				degs.insert(pE2->first + pBi->first);
		}
		std::map<int, int> s_max_map;
		for (auto& deg : degs)
			if (s_max_map[deg.v + 2 * deg.s] < deg.s)
				s_max_map[deg.v + 2 * deg.s] = deg.s;
		
		for (auto [v1, s_max] : s_max_map) {
			Mon1d basis_E2t_sm1, basis_E2t_s, basis_E2t_sp1;
			Poly1d mon_diffs_sp1, mon_diffs_sp2;
			array2d boundaries_s, boundaries_sp1;
			array indices_s, indices_sp1;

			/* s=-1 */
			get_basis_E2t_1(basis_E2, basis_bi, &basis_E2t_sp1, &mon_diffs_sp2, { 0, t, v1 });
			indices_sp1 = range((int)basis_E2t_sp1.size());
			std::sort(indices_sp1.begin(), indices_sp1.end(), [&basis_E2t_sp1](int i1, int i2) {return basis_E2t_sp1[i1] < basis_E2t_sp1[i2]; });
			std::sort(basis_E2t_sp1.begin(), basis_E2t_sp1.end());

			for (int s = 0; s <= s_max; ++s) {
				const Deg& deg = { s, t, v1 - r * s };

				/* Compute cycles and boundaries */
				std::swap(basis_E2t_sm1, basis_E2t_s); std::swap(basis_E2t_s, basis_E2t_sp1);
				std::swap(mon_diffs_sp1, mon_diffs_sp2);
				std::swap(boundaries_s, boundaries_sp1);
				std::swap(indices_s, indices_sp1);
				basis_E2t_sp1.clear();
				mon_diffs_sp2.clear();
				boundaries_sp1.clear();
				get_basis_E2t_1(basis_E2, basis_bi, &basis_E2t_sp1, &mon_diffs_sp2, deg + Deg{ 1, 0, -r });

				indices_sp1 = range((int)basis_E2t_sp1.size());
				std::sort(indices_sp1.begin(), indices_sp1.end(), [&basis_E2t_sp1](int i1, int i2) {return basis_E2t_sp1[i1] < basis_E2t_sp1[i2]; });
				std::sort(basis_E2t_sp1.begin(), basis_E2t_sp1.end());

				array range_ = range((int)basis_E2t_s.size());
				array lead_boundaries;
				for (const auto& base : boundaries_s)
					lead_boundaries.push_back(base.front());
				std::sort(lead_boundaries.begin(), lead_boundaries.end());
				array x = add_vectors(range_, lead_boundaries);
				array2d map_diff;
				for (int xi : x)
					map_diff.push_back(Poly_to_indices(reduce(mon_diffs_sp1[indices_s[xi]], gb), basis_E2t_sp1));

				array2d cycles_s, g_diff;
				set_linear_map(x, map_diff, boundaries_sp1, cycles_s, g_diff);

				auto pNew = basis_H_new.find(deg);
				if (pNew != basis_H_new.end()) {
					/* Compute gb */
					array indices_H = range(int(pNew->second.size()));
					std::sort(indices_H.begin(), indices_H.end(), [&pNew](int a, int b) {return pNew->second[a] < pNew->second[b]; });
					Mon1d mons_new;
					array2d reprs_new;
					for (size_t j = 0; j < indices_H.size(); j++) {
						mons_new.push_back(std::move(pNew->second[indices_H[j]]));
						reprs_new.push_back(residue(boundaries_s, Poly_to_indices(mon_reprs_H_new[deg][indices_H[j]], basis_E2t_s)));
					}

					array2d image, kernel, g, quotient;
					set_linear_map(reprs_new, image, kernel, g);
					array lead_kernel;
					for (array& k : kernel) {
						gb_H.push_back(indices_to_Poly(k, mons_new));
						int index = gb_H.back().front().front().gen;
						if (size_t(index) >= leadings.size())
							leadings.resize(size_t(index) + 1);
						leadings[index].push_back(gb_H.back().front());
						lead_kernel.push_back(k[0]);
					}
					std::sort(lead_kernel.begin(), lead_kernel.end());

					/* Find new generators and add to basis_H */
					quotient = quotient_space(cycles_s, image);
					simplify_space(quotient);
					for (auto& x : quotient) {
						gen_degs_H.push_back(deg);
						reprs_H.push_back(indices_to_Poly(x, basis_E2t_s));
						basis_H[deg].push_back({ {int(gen_degs_H.size()) - 1, 1} });
						mon_reprs_H[deg].push_back(reprs_H.back());
					}

					/* Add to basis_H */
					array index_basis = add_vectors(range(int(mons_new.size())), lead_kernel);
					for (int i : index_basis) {
						basis_H[deg].push_back(std::move(mons_new[i]));
						mon_reprs_H[deg].push_back(mon_reprs_H_new[deg][indices_H[i]]);
					}
				}
				else {
					/* Find new generators and add to basis_H */
					for (array& x : simplify_space(cycles_s)) {
						gen_degs_H.push_back(deg);
						reprs_H.push_back(indices_to_Poly(x, basis_E2t_s));
						basis_H[deg].push_back({ {int(gen_degs_H.size()) - 1, 1} });
						mon_reprs_H[deg].push_back(reprs_H.back());
					}
				}
			}
		}
		for (auto& [deg, basis_H_new_d] : basis_H_new) {
			if (deg.s > s_max_map[deg.v + r * deg.s]) {
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
							const Mon& mon1 = p->second[j];
							const Poly& repr1 = mon_reprs_H[deg1][j];
							if (mon1.empty() || gen_id <= mon1[0].gen) {
								Mon mon(mul(mon1, { {gen_id, 1} }));
								if (gen_id >= (int)leadings.size() || std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(),
									[&mon](Mon _m) { return divides(_m, mon); })) {
									/* Compute the represeting cycle of the base monomial */
									Poly repr = reduce(mul(repr1, reprs_H[gen_id]), gb);
									basis_H_new[deg].push_back(std::move(mon));
									mon_reprs_H_new[deg].push_back(std::move(repr));
								}
							}
						}
					}
				}
			}
		}
	}


	/* Save generators */
	db.save_generators(table_H_prefix + "_generators", gen_degs_H, reprs_H);

	/* Save gb */
	db.save_gb(table_H_prefix + "_relations", gb_H, gen_degs_H);

	/* Save basis */
	db.save_basis(table_H_prefix + "_basis", basis_H, mon_reprs_H);
}

int main_generate_E4t_1(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	db.execute_cmd("DELETE FROM E4t_generators;");
	db.execute_cmd("DELETE FROM E4t_relations;");
	db.execute_cmd("DELETE FROM E4t_basis;");

	auto start = std::chrono::system_clock::now();

	int t_max = 74;
	generate_E4t(db, t_max);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Elapsed time: " << elapsed.count() << "s\n";

	return 0;
}