#include "main_E4t.h"
#include "sqlite3/sqlite3.h"
#include "database.h"
#include "linalg.h"

/* assume that the sequence map_gen_id is increasing */
Poly reindex(const Poly& poly, const array& map_gen_id)
{
	Poly result;
	for (const Mon& m : poly) {
		Mon m1;
		for (GenPow ge : m)
			m1.push_back({ map_gen_id[ge.gen], ge.exp });
		result.push_back(m1);
	}
	return result;
}

std::map<Deg, DgaBasis1> load_dga_basis(const Database& db, const std::string& table_name, int r, int t_max)
{
	std::map<Deg, DgaBasis1> basis;
	Statement stmt(db, "SELECT mon, diff, s, t, v FROM " + table_name + " WHERE t<=" + std::to_string(t_max) + " ORDER BY mon_id;");
	int prev_t = 0;
	std::map<Deg, array2d> mon_diffs;
	while (stmt.step() == SQLITE_ROW) {
		const Deg d = { stmt.column_int(2), stmt.column_int(3), stmt.column_int(4) };
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
			Deg d1 = d + Deg{ 1, 0, -r };
			basis_d.diffs[i] = basis.find(d1) != basis.end() ? indices_to_Poly(mon_diffs[d][i], basis.at(d1).basis) : Poly{};
		}
	}
	std::cout << "dga_basis loaded from " << table_name << ", # of degs=" << basis.size() << '\n';
	return basis;
}

std::map<Deg, DgaBasis1> load_basis_X(const Database& db, const std::string& table_name, int t_max, int r)
{
	std::map<Deg, DgaBasis1> basis;
	Statement stmt(db, "SELECT mon, diff, s, t, v FROM " + table_name + " WHERE t<=" + std::to_string(t_max) + " ORDER BY mon_id;");
	int prev_t = 0;
	while (stmt.step() == SQLITE_ROW) {
		Deg d = { stmt.column_int(2), stmt.column_int(3), stmt.column_int(4) };
		basis[d].basis.push_back(str_to_Mon(stmt.column_str(0)));
		basis[d].diffs.push_back(str_to_Poly(stmt.column_str(1)));
	}
	std::cout << "basis loaded from " << table_name << ", size=" << basis.size() << '\n';
	return basis;
}

void load_y(const Database& db, const std::string& table_name, Poly1d& y, array& t_y)
{
	Statement stmt(db, "SELECT y, t FROM " + table_name + " ORDER BY t;");
	while (stmt.step() == SQLITE_ROW) {
		y.push_back(str_to_Poly(stmt.column_str(0)));
		t_y.push_back(stmt.column_int(1));
	}
	std::cout << "y loaded from " << table_name << ", size=" << y.size() << '\n';
}

void save_y(const Database& db, const std::string& table_name, const Poly2d& y_t, int t)
{
	Statement stmt_update_relations(db, "INSERT INTO " + table_name + " (y, t) VALUES (?1, ?2);");

	for (size_t i = 0; i < y_t.size(); ++i) {
		stmt_update_relations.bind_str(1, Poly_to_str(y_t[i][0]));
		stmt_update_relations.bind_int(2, t);
		stmt_update_relations.step_and_reset();
	}
#ifdef DATABASE_SAVE_LOGGING
	std::cout << y_t.size() << " y's are inserted!\n";
#endif
}

void save_map_gen_id(const Database& db, const std::string& table_name, const array& map_gen_id, int i_start)
{
	Statement stmt_update_relations(db, "INSERT INTO " + table_name + " (gen_id) VALUES (?1);");

	for (size_t i = i_start; i < map_gen_id.size(); ++i) {
		stmt_update_relations.bind_int(1, map_gen_id[i]);
		stmt_update_relations.step_and_reset();
	}
#ifdef DATABASE_SAVE_LOGGING
	std::cout << map_gen_id.size() - i_start << " new gen_id's are inserted!\n";
#endif
}

void SaveGb(const Database& db, const std::string& table_name, const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t)
{
	Statement stmt_update_relations(db, "INSERT INTO " + table_name + " (leading_term, basis, t) VALUES (?1, ?2, ?3);"); ////

	for (int i = (int)gb.size() - 1; i >= 0; --i) {
		int t1 = get_deg(gb[i].front(), gen_degs, gen_degs1);
		if (t1 != t)
			break;
		stmt_update_relations.bind_str(1, Mon_to_str(gb[i].front()));
		stmt_update_relations.bind_str(2, Poly_to_str(gb[i].begin() + 1, gb[i].end()));
		stmt_update_relations.bind_int(3, t1); ////
		stmt_update_relations.step_and_reset();
	}
#ifdef DATABASE_SAVE_LOGGING
	std::cout << gb.size() << " relations are inserted into " + table_name + "!\n";
#endif
}

void get_basis_B(const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X, Mon1d& basis_B, const Deg& deg)
{
	for (auto pX = basis_X.rbegin(); pX != basis_X.rend(); ++pX) {
		const Deg& d2 = pX->first;
		if (d2.s <= deg.s && d2.t <= deg.t && d2.v <= deg.v) {
			const Deg d1 = deg - d2;
			auto pA = basis_A.find(d1);
			if (pA != basis_A.end())
				for (size_t i = 0; i < pA->second.basis.size(); ++i)
					for (size_t j = 0; j < pX->second.basis.size(); ++j)
						basis_B.push_back(mul(pA->second.basis[i], pX->second.basis[j]));
		}

	}
}

void get_basis_with_diff_B(const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X, Mon1d& basis_B, Poly1d& mon_diffs_B, const Deg& deg)
{
	for (auto pX = basis_X.rbegin(); pX != basis_X.rend(); ++pX) {
		const Deg& d2 = pX->first;
		if (d2.s <= deg.s && d2.t <= deg.t && d2.v <= deg.v) {
			const Deg d1 = deg - d2;
			auto pA = basis_A.find(d1);
			if (pA != basis_A.end())
				for (size_t i = 0; i < pA->second.basis.size(); ++i)
					for (size_t j = 0; j < pX->second.basis.size(); ++j) {
						basis_B.push_back(mul(pA->second.basis[i], pX->second.basis[j]));
						Poly diff = add(mul(pA->second.diffs[i], pX->second.basis[j]), mul(pX->second.diffs[j], pA->second.basis[i]));
						mon_diffs_B.push_back(std::move(diff));
					}
		}

	}
}

Poly d_inv(const Poly& poly, const std::vector<Deg>& gen_degs, const Poly1d& diffs, const grbn::GbWithCache& gb, const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X)
{
	if (poly.empty())
		return {};
	Deg deg_poly = get_deg(poly[0], gen_degs);
	Deg deg_result = deg_poly - Deg{ 1, 0, -2 };
	Mon1d basis_in_poly;
	Mon1d basis_in_result;
	get_basis_B(basis_A, basis_X, basis_in_poly, deg_poly);
	get_basis_B(basis_A, basis_X, basis_in_result, deg_result);
	std::sort(basis_in_poly.begin(), basis_in_poly.end());
	std::sort(basis_in_result.begin(), basis_in_result.end());
	array2d map_diff;
	for (const Mon& mon : basis_in_result)
		map_diff.push_back(Poly_to_indices(grbn::Reduce(get_diff(mon, diffs), gb), basis_in_poly));
	array2d image, kernel, g;
	lina::SetLinearMap(map_diff, image, kernel, g);
	return indices_to_Poly(lina::GetImage(image, g, Poly_to_indices(poly, basis_in_poly)), basis_in_result);
}

Poly proj(const Poly& poly, const std::vector<Deg>& gen_degs, const Poly1d& gen_diffs, const grbn::GbWithCache& gb, const std::map<Deg, DgaBasis1>& basis_A,
	const std::map<Deg, DgaBasis1>& basis_X, const Poly1d& gen_reprs_H, std::map<Deg, Mon1d>& basis_H) //
{
	if (poly.empty())
		return {};
	Deg d = get_deg(poly[0], gen_degs);
	Deg d1 = d - Deg{ 1, 0, -2 };
	Mon1d basis_d;
	Mon1d basis_d1;
	get_basis_B(basis_A, basis_X, basis_d, d);
	get_basis_B(basis_A, basis_X, basis_d1, d1);
	std::sort(basis_d.begin(), basis_d.end());
	std::sort(basis_d1.begin(), basis_d1.end());

	array2d map_diff;
	for (const Mon& mon : basis_d1)
		map_diff.push_back(Poly_to_indices(grbn::Reduce(get_diff(mon, gen_diffs), gb), basis_d)); //
	array2d image = lina::GetSpace(map_diff);

	array2d map_repr;
	for (const Mon& mon : basis_H[d])
		map_repr.push_back(lina::Residue(image, Poly_to_indices(get_image({ mon }, gen_reprs_H, gb), basis_d)));
	array2d image1, kernel1, g1;
	lina::SetLinearMap(map_repr, image1, kernel1, g1);

	array poly_ind_mod_boundary = lina::Residue(image, Poly_to_indices(poly, basis_d));
	return indices_to_Poly(lina::GetImage(image1, g1, std::move(poly_ind_mod_boundary)), basis_H[d]);
}