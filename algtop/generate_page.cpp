#include "main.h"

/* Generate the homology of the E_r page */
void generate_next_page(const Database& db, const std::string& table_prefix, const std::string& table_H_prefix, int r)
{
	array3d gb = db.load_gb(table_prefix + "_relations");
	std::map<Deg, array2d> basis = db.load_basis(table_prefix + "_basis");
	std::map<Deg, BasisComplex> basis_ss = db.load_basis_ss(table_prefix + "_ss", r);

	std::vector<Deg> gen_degs_H;
	array2d reprs_H;
	array3d gb_H;
	
	array3d leadings;
	std::map<Deg, array2d> basis_H;
	basis[Deg{ 0, 0, 0 }].push_back({});
	std::map<Deg, array2d> mon_reprs_H;
	basis[Deg{ 0, 0, 0 }].push_back({});
	std::map<Deg, array2d> basis_H_new;
	std::map<Deg, array2d> mon_reprs_H_new;

	int t_max = basis_ss.crbegin()->first.t;
	for (int t = 1; t <= t_max; ++t){
		std::cout << "generate next page, t=" << t << "          \r";
		auto p1_ss = basis_ss.lower_bound(Deg{ 0, t, 0 });
		auto p2_ss = basis_ss.lower_bound(Deg{ 0, t + 1, 0 });
		for (auto p_ss = p1_ss; p_ss != p2_ss; ++p_ss) {
			const Deg& deg = p_ss->first;
			auto pNew = basis_H_new.find(deg);
			if (pNew != basis_H_new.end()) {
				/* Compute gb */
				array indices = range(int(pNew->second.size()));
				std::sort(indices.begin(), indices.end(), [&pNew](int a, int b) {return cmp_mons(pNew->second[a], pNew->second[b]); });
				array2d mons_new, reprs_new;
				for (size_t j = 0; j < indices.size(); j++) {
					mons_new.push_back(std::move(pNew->second[indices[j]]));
					reprs_new.push_back(std::move(mon_reprs_H_new[deg][indices[j]]));
				}

				array2d image, kernel, g, quotient;
				set_linear_map(reprs_new, image, kernel, g);
				array lead_kernel;
				for (array& p_indices : kernel) {
					gb_H.push_back(indices_to_poly(p_indices, mons_new));
					int index = gb_H.back().front().front();
					if (size_t(index) >= leadings.size())
						leadings.resize(size_t(index) + 1);
					leadings[index].push_back(gb_H.back().front());
					lead_kernel.push_back(p_indices[0]);
				}
				std::sort(lead_kernel.begin(), lead_kernel.end());

				/* Find new generators and add to basis_H */
				quotient = quotient_space(p_ss->second.kernel, image);
				simplify_space(quotient);
				for (auto& x : quotient) {
					gen_degs_H.push_back(deg);
					reprs_H.push_back(x);
					basis_H[deg].push_back({ int(gen_degs_H.size()) - 1, 1 });
					mon_reprs_H[deg].push_back(std::move(x));
				}

				/* Add to basis_H */
				array index_basis = add_vectors(range(int(mons_new.size())), lead_kernel);
				for (int i : index_basis) {
					basis_H[deg].push_back(std::move(mons_new[i]));
					mon_reprs_H[deg].push_back(std::move(reprs_new[i]));
				}
			}
			else {
				/* Find new generators and add to basis_H */
				for (array& x : simplify_space(p_ss->second.kernel)) {
					gen_degs_H.push_back(deg);
					reprs_H.push_back(x);
					basis_H[deg].push_back({ int(gen_degs_H.size()) - 1, 1 });
					mon_reprs_H[deg].push_back(std::move(x));
				}
			}
		}
		for (auto& [deg, basis_H_new_d] : basis_H_new) {
			if (basis_ss.find(deg) == basis_ss.end()) {
				/* Add to gb_H */
				for (auto& mon : basis_H_new_d) {
					gb_H.push_back({{std::move(mon)}});
					int index = gb_H.back().front().front();
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
							if (mon1.empty() || gen_id <= mon1[0]) {
								array mon(mul(mon1, { gen_id, 1 }));
								if (gen_id >= (int)leadings.size() || std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(),
									[mon](array _m) { return divides(_m, mon); })) {
									/* Compute the represeting cycle of the base monomial */
									auto poly1 = indices_to_poly(repr1, basis[deg1]);
									auto gen_repr = indices_to_poly(reprs_H[gen_id], basis[g_deg]);
									array2d repr = reduce(mul(poly1, gen_repr), gb);
									array repr_indices = residue(basis_ss[deg].boundary, poly_to_indices(repr, basis[deg]));
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

	/* Save generators */
	array3d reprs_poly_H;
	for (size_t i = 0; i < (int)reprs_H.size(); ++i)
		reprs_poly_H.push_back(indices_to_poly(reprs_H[i], basis[gen_degs_H[i]]));
	db.save_generators(table_H_prefix + "_generators", gen_degs_H, reprs_poly_H);

	/* Save gb */
	db.save_gb(table_H_prefix + "_relations", gb_H, gen_degs_H);

	/* Save basis */
	db.save_basis(table_H_prefix + "_basis", basis_H, mon_reprs_H);
}

int main_generate_next_page(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	std::string table_prefix = "E2", table_H_prefix = "E4";
	generate_next_page(db, table_prefix, table_H_prefix, 2);
	return 0;
}