#include "database.h"
#include <sstream>

struct Basis_d
{
	Deg deg;
	array2d mons;
};

struct BasisSS_d
{
	Deg deg;
	array2d boundary;
	array2d kernel;
};

struct Generator_H
{
	Deg deg;
	array repr;

	Generator_H(const Deg& deg_, array&& repr_) : deg(deg_), repr(repr_) {};
};

struct Basis_H_d
{
	Deg deg;
	array2d mons;
	array2d reprs;
};

inline array2d indices_to_poly(const array& indices, const std::vector<Basis_d>& basis, const Deg& deg)
{
	array2d result;
	const auto& basis_d = get_item_by_deg(basis, deg);
	for (int i : indices)
		result.push_back(basis_d.mons[i]);
	return result;
}

inline array2d indices_to_poly(const array& indices, const std::vector<array>& basis)
{
	array2d result;
	for (int i : indices)
		result.push_back(basis[i]);
	return result;
}

inline int get_index(const std::vector<array>& basis, const array& mon)
{
	auto index = std::lower_bound(basis.begin(), basis.end(), mon, cmp_mons);
	return int(index - basis.begin());
}

inline array poly_to_indices(const array2d& poly, const std::vector<Basis_d>& basis, const Deg& deg)
{
	array result;
	const auto& basis_d = get_item_by_deg(basis, deg);
	for (const array& mon : poly)
		result.push_back(get_index(basis_d.mons, mon));
	return result;
}

void load_basis(sqlite3* conn, const std::string& table_name_basis, std::vector<Basis_d>& basis)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT mon, s, t, v FROM ") + table_name_basis + " ORDER BY mon_id;";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	Deg prev_deg = { -1, -1, -1 };
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		Deg deg = { sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3) };
		if (deg == prev_deg)
			basis.back().mons.push_back(str_to_array(sqlite3_column_str(stmt, 0)));
		else
			basis.push_back(Basis_d({ deg, {str_to_array(sqlite3_column_str(stmt, 0))} }));
		prev_deg = deg;
	}
	sqlite3_finalize(stmt);
}

void load_basis_ss(sqlite3* conn, const std::string& table_name_ss, std::vector<BasisSS_d>& basis_ss, int r)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT base, level, s, t, v FROM ") + table_name_ss + " ;";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	Deg prev_deg = { -1, -1, -1 };
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		Deg deg = { sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3), sqlite3_column_int(stmt, 4) };
		int level = sqlite3_column_int(stmt, 1);
		if (deg == prev_deg) {
			if (level <= r)
				basis_ss.back().boundary.push_back(str_to_array(sqlite3_column_str(stmt, 0)));
			else if (level <= N - r)
				basis_ss.back().kernel.push_back(str_to_array(sqlite3_column_str(stmt, 0)));
		}
		else {
			if (level <= r)
				basis_ss.push_back(BasisSS_d({ deg, {str_to_array(sqlite3_column_str(stmt, 0))}, {} }));
			else if (level <= N - r)
				basis_ss.push_back(BasisSS_d({ deg, {}, {str_to_array(sqlite3_column_str(stmt, 0))} }));
		}
		prev_deg = deg;
	}
	sqlite3_finalize(stmt);
}

/* Generate the homology of the E_r page */
void generate_next_page(sqlite3* conn, const std::string& table_prefix, const std::string& table_H_prefix, int r)
{
	/* Load relations */
	array3d relations;
	load_relations(conn, table_prefix + "_relations", relations);
	std::cout << "relations loaded! Size=" << relations.size() << '\n';

	/* Load basis */
	std::vector<Basis_d> basis;
	load_basis(conn, table_prefix + "_basis", basis);
	std::cout << "basis loaded! Size=" << basis.size() << '\n';

	/* Load basis_ss */
	std::vector<BasisSS_d> basis_ss;
	load_basis_ss(conn, table_prefix + "_ss", basis_ss, r);
	std::cout << "basis_ss loaded! Size=" << basis_ss.size() << '\n';

	std::vector<Generator_H> generators_H;
	array3d relations_H;
	array3d leadings;
	std::vector<Basis_H_d> basis_H = { {{0, 0, 0}, {{}}, {{}}} };
	std::vector<size_t> num_degs_t = { 1 };
	std::vector<Basis_H_d> basis_H_new;

	int t_max = basis.back().deg.t;
	int i = 1;
	for (int t = 1; t <= t_max; t++) {
		std::cout << "t=" << t << '\n';
		for (; i < basis_ss.size() && basis_ss[i].deg.t == t; i++) {
			auto& basis_ss_d = basis_ss[i];
			const auto& deg = basis_ss_d.deg;
			basis_H.push_back(Basis_H_d({ deg, {}, {} }));
			if (deg_in(deg, basis_H_new)) {
				auto& basis_H_new_d = get_item_by_deg(basis_H_new, deg);
				array2d image, kernel, quotient;

				/* Compute relations */
				get_image_kernel(basis_H_new_d.reprs, image, kernel);
				array lead_kernel;
				for (array& m : kernel) {
					std::sort(m.begin(), m.end(), [&basis_H_new_d](int a, int b) {return cmp_mons(basis_H_new_d.mons[a], basis_H_new_d.mons[b]); });
					array2d rel = indices_to_poly(m, basis_H_new_d.mons);
					relations_H.push_back(std::move(rel));
					int index = relations_H.back()[0][0];
					if (size_t(index) >= leadings.size())
						leadings.resize(size_t(index) + 1);
					leadings[index].push_back(relations_H.back()[0]);
					lead_kernel.push_back(m[0]);
				}
				std::sort(lead_kernel.begin(), lead_kernel.end());

				/* Find new generators and add to basis_H */
				quotient = quotient_space(basis_ss_d.kernel, image);
				simplify_space(quotient);
				for (auto& x : quotient) {
					generators_H.emplace_back(deg, array(x));
					basis_H.back().mons.push_back({ int(generators_H.size()) - 1, 1 });
					basis_H.back().reprs.push_back(std::move(x));
				}

				/* Add to basis_H */
				array index_basis = add_vectors(range(int(basis_H_new_d.mons.size())), lead_kernel);
				for (int i : index_basis) {
					basis_H.back().mons.push_back(std::move(basis_H_new_d.mons[i]));
					basis_H.back().reprs.push_back(std::move(basis_H_new_d.reprs[i]));
				}
			}
			else {
				/* Find new generators and add to basis_H */
				for (array& x : simplify_space(basis_ss_d.kernel)) {
					generators_H.emplace_back(deg, array(x));
					basis_H.back().mons.push_back({ int(generators_H.size()) - 1, 1 });
					basis_H.back().reprs.push_back(std::move(x));
				}
			}
		}
		for (auto& basis_H_new_d : basis_H_new) {
			if (!deg_in(basis_H_new_d.deg, basis_ss)) {
				/* Add to relations_H */
				for (auto& mon : basis_H_new_d.mons) {
					relations_H.push_back({ {std::move(mon)} });
					int index = relations_H.back()[0][0];
					if (size_t(index) >= leadings.size())
						leadings.resize(size_t(index) + 1);
					leadings[index].push_back(relations_H.back()[0]);
				}
			}
		}
		num_degs_t.push_back(basis_H.size());

		/* Compute new basis of degree t + 1 */
		if (t < t_max) {
			basis_H_new.clear();
			for (int id = 0; size_t(id) < generators_H.size(); id++) {
				const auto& gen = generators_H[id];
				int t1 = t + 1 - gen.deg.t;
				if (t1 >= 0) {
					size_t index1 = t1 > 0 ? num_degs_t[size_t(t1) - 1] : 0;
					size_t index2 = num_degs_t[t1];
					for (size_t i = index1; i < index2; i++) {
						Deg deg1 = basis_H[i].deg;
						Deg deg = deg1 + gen.deg;
						for (size_t j = 0; j < basis_H[i].mons.size(); j++) {
							const auto& mon1 = basis_H[i].mons[j];
							const auto& repr1 = basis_H[i].reprs[j];
							if (mon1.empty() || id <= mon1[0]) {
								array mon(mul(mon1, { id, 1 }));
								if ((size_t)id >= leadings.size() || std::none_of(leadings[id].begin(), leadings[id].end(),
									[mon](array _m) { return divides(_m, mon); })) {
									/* Compute the represeting cycle of the base monomial */
									auto p1 = indices_to_poly(repr1, basis, deg1);
									auto p2 = indices_to_poly(gen.repr, basis, gen.deg);
									array2d prod = mul(p1, p2);
									prod = reduce(prod, relations);
									array prod_indices = poly_to_indices(prod, basis, deg);
									const auto& basis_ss_d = get_item_by_deg(basis_ss, deg);
									array prod_indices1 = residue(basis_ss_d.boundary, std::move(prod_indices));
									auto pEle = std::lower_bound(basis_H_new.begin(), basis_H_new.end(), deg, [](const Basis_H_d& ele, const Deg& d) {return ele.deg < d; });
									if (pEle == basis_H_new.end() || pEle->deg != deg)
										pEle = basis_H_new.insert(pEle, Basis_H_d({ deg, {}, {} }));
									pEle->mons.push_back(std::move(mon));
									pEle->reprs.push_back(std::move(prod_indices1));
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
	for (int i = 0; i < int(generators_H.size()); i++) {
		sqlite3_bind_int(stmt_update_generators, 1, i);
		sqlite3_bind_str(stmt_update_generators, 2, array_to_str(generators_H[i].repr));
		sqlite3_bind_int(stmt_update_generators, 3, generators_H[i].deg.s);
		sqlite3_bind_int(stmt_update_generators, 4, generators_H[i].deg.t);
		sqlite3_bind_int(stmt_update_generators, 5, generators_H[i].deg.v);
		sqlite3_step(stmt_update_generators);
		sqlite3_reset(stmt_update_generators);
	}
	execute_cmd(conn, "END TRANSACTION");
	std::cout << generators_H.size() << " Generators are inserted!\n";

	sqlite3_finalize(stmt_update_generators);

	/* Save relations */
	sqlite3_stmt* stmt_update_relations;
	std::string cmd_update_relations = std::string("INSERT INTO ") + table_H_prefix + "_relations (leading_term, basis) VALUES (?1, ?2);"; // Insert s, t, v
	sqlite3_prepare_v100(conn, cmd_update_relations, &stmt_update_relations);

	execute_cmd(conn, "BEGIN TRANSACTION");
	for (int i = 0; i < int(relations_H.size()); i++) {
		sqlite3_bind_str(stmt_update_relations, 1, array_to_str(relations_H[i][0]));
		sqlite3_bind_str(stmt_update_relations, 2, array2d_to_str(relations_H[i].begin() + 1, relations_H[i].end()));
		sqlite3_step(stmt_update_relations);
		sqlite3_reset(stmt_update_relations);
	}
	execute_cmd(conn, "END TRANSACTION");
	std::cout << relations_H.size() << " relations are inserted!\n";

	sqlite3_finalize(stmt_update_relations);

	/* Save basis */
	sqlite3_stmt* stmt_update_basis;
	std::string cmd_update_basis = std::string("INSERT INTO ") + table_H_prefix + "_basis (mon_id, mon, repr, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5, ?6);"; // Insert s, t, v
	sqlite3_prepare_v100(conn, cmd_update_basis, &stmt_update_basis);

	execute_cmd(conn, "BEGIN TRANSACTION");
	size_t mon_id = 0;
	for (size_t i = 0; i < basis_H.size(); i++) {
		array indices = range(int(basis_H[i].mons.size()));
		std::sort(indices.begin(), indices.end(), [&basis_H, i](int a, int b) {return cmp_mons(basis_H[i].mons[a], basis_H[i].mons[b]); });
		for (size_t j = 0; j < indices.size(); j++) {
			sqlite3_bind_int(stmt_update_basis, 1, mon_id);
			sqlite3_bind_str(stmt_update_basis, 2, array_to_str(basis_H[i].mons[indices[j]]));
			sqlite3_bind_str(stmt_update_basis, 3, array_to_str(basis_H[i].reprs[indices[j]]));
			sqlite3_bind_int(stmt_update_basis, 4, basis_H[i].deg.s);
			sqlite3_bind_int(stmt_update_basis, 5, basis_H[i].deg.t);
			sqlite3_bind_int(stmt_update_basis, 6, basis_H[i].deg.v);
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
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\ss_test.db)", &conn);

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
	return 0;
}