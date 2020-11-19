#include "main.h"

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

void load_y(const Database& db, const std::string& table_name, Poly1d& y, array& t_y)
{
	Statement stmt;
	stmt.init(db, "SELECT y, t FROM " + table_name + " ORDER BY t;");
	while (stmt.step() == SQLITE_ROW) {
		y.push_back(str_to_Poly(stmt.column_str(0)));
		t_y.push_back(stmt.column_int(1));
	}
	std::cout << "y loaded from " << table_name << ", size=" << y.size() << '\n';
}

void save_y(const Database& db, const std::string& table_name, const Poly2d& y_t, int t)
{
	Statement stmt_update_relations;
	stmt_update_relations.init(db, "INSERT INTO " + table_name + " (y, t) VALUES (?1, ?2);");

	for (size_t i = 0; i < y_t.size(); ++i) {
		stmt_update_relations.bind_str(1, Poly_to_str(y_t[i][0]));
		stmt_update_relations.bind_int(2, t);
		stmt_update_relations.step_and_reset();
	}
	std::cout << y_t.size() << " y's are inserted!\n";
}

void SaveGb(const Database& db, const std::string& table_name, const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t)
{
	Statement stmt_update_relations;
	stmt_update_relations.init(db, "INSERT INTO " + table_name + " (leading_term, basis) VALUES (?1, ?2);");

	for (int i = (int)gb.size() - 1; i >= 0; --i) {
		int t1 = get_deg(gb[i].front(), gen_degs, gen_degs1);
		if (t1 != t)
			break;
		stmt_update_relations.bind_str(1, Mon_to_str(gb[i].front()));
		stmt_update_relations.bind_str(2, Poly_to_str(gb[i].begin() + 1, gb[i].end()));
		stmt_update_relations.step_and_reset();
	}
	std::cout << gb.size() << " relations are inserted into " + table_name + "!\n";
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

Poly d_inv(const Poly& poly, const std::vector<Deg>& gen_degs, const Poly1d& diffs, const Poly1d& gb, const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X)
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
		map_diff.push_back(Poly_to_indices(reduce(get_diff(mon, diffs), gb), basis_in_poly));
	array2d image, kernel, g;
	set_linear_map(map_diff, image, kernel, g);
	return indices_to_Poly(get_image(image, g, Poly_to_indices(poly, basis_in_poly)), basis_in_result);
}

Poly proj(const Poly& poly, const std::vector<Deg>& gen_degs, const Poly1d& gen_diffs, const Poly1d& gb, const std::map<Deg, DgaBasis1>& basis_A,
	const std::map<Deg, DgaBasis1>& basis_X, const Poly1d& gen_reprs, std::map<Deg, Mon1d>& basis_H)
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
		map_diff.push_back(Poly_to_indices(reduce(get_diff(mon, gen_diffs), gb), basis_in_poly));
	array2d image, kernel, g;
	set_linear_map(map_diff, image, kernel, g);

	array2d map_repr;
	for (const Mon& mon : basis_H[deg_poly])
		map_repr.push_back(residue(image, Poly_to_indices(get_repr({ mon }, gen_reprs, gb), basis_in_poly)));
	array2d image1, kernel1, g1;
	set_linear_map(map_repr, image1, kernel1, g1);

	return indices_to_Poly(get_image(image1, g1, Poly_to_indices(poly, basis_in_poly)), basis_H[deg_poly]);
}