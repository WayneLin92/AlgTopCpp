#include "database.h"

/********** STRUCTS AND CLASSES **********/

struct BasisSS1
{
	array2d boundary;
	array2d kernel;
};

/********** FUNCTIONS **********/

array2d indices_to_poly(const array& indices, const array2d& basis)
{
	array2d result;
	for (int i : indices)
		result.push_back(basis[i]);
	return result;
}

array poly_to_indices(const array2d& poly, const array2d& basis)
{
	array result;
	for (const array& mon : poly)
		result.push_back(get_index(basis, mon));
	return result;
}

void load_basis_ss(sqlite3* conn, const std::string& table_name_ss, std::map<Deg, BasisSS1>& basis_ss, int r)
{
	sqlite3_stmt* stmt;
	std::string cmd = "SELECT base, level, s, t, v FROM " + table_name_ss + " ;";
	sqlite3_prepare_v100(conn, cmd, &stmt);

	while (sqlite3_step(stmt) == SQLITE_ROW) {
		Deg deg = { sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4) };
		int level = sqlite3_column_int(stmt, 1);
		if (level <= r)
			basis_ss[deg].boundary.push_back(str_to_array(sqlite3_column_str(stmt, 0)));
		else if (level <= T_MAX - r)
			basis_ss[deg].kernel.push_back(str_to_array(sqlite3_column_str(stmt, 0)));
	}
	sqlite3_finalize(stmt);
}

/* Generate the homology of the E_r page */
void generate_next_page(sqlite3* conn, const std::string& table_prefix, const std::string& table_H_prefix, int r)
{
	/* Load gb */
	array3d gb;
	load_gb(conn, table_prefix + "_relations", gb);
	std::cout << "gb loaded! Size=" << gb.size() << '\n';

	/* Load basis */
	std::map<Deg, array2d> basis;
	load_basis(conn, table_prefix + "_basis", basis);
	std::cout << "basis loaded! Size=" << basis.size() << '\n';

	/* Load basis_ss */
	std::map<Deg, BasisSS1> basis_ss;
	load_basis_ss(conn, table_prefix + "_ss", basis_ss, r);
	std::cout << "basis_ss loaded! Size=" << basis_ss.size() << '\n';

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
	int prev_t = 0;
	for (int t = 1; t <= t_max; ++t){
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
	sqlite3_stmt* stmt_update_generators;
	std::string cmd_update_generators = std::string("INSERT INTO ") + table_H_prefix + "_generators (gen_id, repr, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5);";
	sqlite3_prepare_v100(conn, cmd_update_generators, &stmt_update_generators);

	execute_cmd(conn, "BEGIN TRANSACTION");
	for (int i = 0; i < (int)gen_degs_H.size(); ++i) {
		sqlite3_bind_int(stmt_update_generators, 1, i);
		sqlite3_bind_str(stmt_update_generators, 2, array2d_to_str(indices_to_poly(reprs_H[i], basis[gen_degs_H[i]])));
		sqlite3_bind_int(stmt_update_generators, 3, gen_degs_H[i].s);
		sqlite3_bind_int(stmt_update_generators, 4, gen_degs_H[i].t);
		sqlite3_bind_int(stmt_update_generators, 5, gen_degs_H[i].v);
		sqlite3_step(stmt_update_generators);
		sqlite3_reset(stmt_update_generators);
	}
	execute_cmd(conn, "END TRANSACTION");
	std::cout << gen_degs_H.size() << " Generators are inserted!\n";

	sqlite3_finalize(stmt_update_generators);

	/* Save gb */
	sqlite3_stmt* stmt_update_relations;
	std::string cmd_update_relations = std::string("INSERT INTO ") + table_H_prefix + "_relations (leading_term, basis, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5);"; // Insert s, t, v
	sqlite3_prepare_v100(conn, cmd_update_relations, &stmt_update_relations);

	execute_cmd(conn, "BEGIN TRANSACTION");
	for (int i = 0; i < int(gb_H.size()); ++i) {
		Deg deg = get_deg(gb_H[i], gen_degs_H);
		sqlite3_bind_str(stmt_update_relations, 1, array_to_str(gb_H[i].front()));
		sqlite3_bind_str(stmt_update_relations, 2, array2d_to_str(gb_H[i].begin() + 1, gb_H[i].end()));
		sqlite3_bind_int(stmt_update_relations, 3, deg.s);
		sqlite3_bind_int(stmt_update_relations, 4, deg.t);
		sqlite3_bind_int(stmt_update_relations, 5, deg.v);
		sqlite3_step(stmt_update_relations);
		sqlite3_reset(stmt_update_relations);
	}
	execute_cmd(conn, "END TRANSACTION");
	std::cout << gb_H.size() << " relations are inserted!\n";

	sqlite3_finalize(stmt_update_relations);

	/* Save basis */
	sqlite3_stmt* stmt_update_basis;
	std::string cmd_update_basis = std::string("INSERT INTO ") + table_H_prefix + "_basis (mon_id, mon, repr, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5, ?6);"; // Insert s, t, v
	sqlite3_prepare_v100(conn, cmd_update_basis, &stmt_update_basis);

	execute_cmd(conn, "BEGIN TRANSACTION");
	int mon_id = 0;
	for (auto p = basis_H.begin(); p != basis_H.end(); ++p) {
		for (size_t i = 0; i < p->second.size(); ++i) {
			const Deg& deg = p->first;
			sqlite3_bind_int(stmt_update_basis, 1, mon_id);
			sqlite3_bind_str(stmt_update_basis, 2, array_to_str(p->second[i]));
			sqlite3_bind_str(stmt_update_basis, 3, array_to_str(mon_reprs_H[deg][i]));
			sqlite3_bind_int(stmt_update_basis, 4, deg.s);
			sqlite3_bind_int(stmt_update_basis, 5, deg.t);
			sqlite3_bind_int(stmt_update_basis, 6, deg.v);
			sqlite3_step(stmt_update_basis);
			sqlite3_reset(stmt_update_basis);
			mon_id++;
		}
	}
	execute_cmd(conn, "END TRANSACTION");
	std::cout << basis_H.size() << " basis are inserted!\n";

	sqlite3_finalize(stmt_update_basis);
}

int main_generate_next_page(int argc, char** argv)
{
	sqlite3* conn;
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)", &conn);

	std::string table_prefix, table_H_prefix;
	if (argc == 1) {
		table_prefix = "E2";
		table_H_prefix = "E4";
	}
	else {
		table_prefix = argv[1];
		table_H_prefix = argv[2];
	}

	generate_next_page(conn, table_prefix, table_H_prefix, 2);

	sqlite3_close(conn);
	return 0;
}