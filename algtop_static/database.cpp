#include "database.h"
#include "sqlite3/sqlite3.h"
#include "myio.h"
#include <iostream>
#include <sstream>

Database::~Database()
{
	sqlite3_close(conn_);
}

void Database::init(const char* filename)
{
	if (sqlite3_open(filename, &conn_) != SQLITE_OK)
		throw "8de81e80";
}

void Database::execute_cmd(const std::string& sql) const
{
	sqlite3_stmt* stmt;
	sqlite3_prepare_v100(sql, &stmt);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
}

int Database::get_int(const std::string& sql) const
{
	Statement stmt;
	stmt.init(*this, sql);
	if (stmt.step() == SQLITE_ROW)
		if (stmt.column_type(0) == SQLITE_INTEGER)
			return stmt.column_int(0);
	throw "effbf28c";
}

array Database::get_ints(const std::string& table_name, const std::string& column_name, const std::string& conditions) const
{
	array result;
	Statement stmt;
	stmt.init(*this, "SELECT " + column_name + " FROM " + table_name + ' ' + conditions + ';');
	while (stmt.step() == SQLITE_ROW)
		result.emplace_back(stmt.column_int(0));
	std::cout << column_name << " loaded from " << table_name << ", size=" << result.size() << '\n';
	return result;
}

void Database::sqlite3_prepare_v100(const char* zSql, sqlite3_stmt** ppStmt) const
{
	if (sqlite3_prepare_v2(conn_, zSql, int(strlen(zSql)) + 1, ppStmt, NULL) != SQLITE_OK)
		throw "bce2dcfe";
}

void Database::sqlite3_prepare_v100(const std::string& sql, sqlite3_stmt** ppStmt) const
{
	if (sqlite3_prepare_v2(conn_, sql.c_str(), int(sql.size()) + 1, ppStmt, NULL) != SQLITE_OK)
		throw "da6ab7f6";
}

std::vector<Deg> Database::load_gen_degs(const std::string& table_name) const
{
	std::vector<Deg> gen_degs;
	Statement stmt;
	stmt.init(*this, "SELECT s, t, v FROM " + table_name + " ORDER BY gen_id;");
	while (stmt.step() == SQLITE_ROW)
		gen_degs.emplace_back(stmt.column_int(0), stmt.column_int(1), stmt.column_int(2));
	std::cout << "gen_degs loaded from " << table_name << ", size=" << gen_degs.size() << '\n';
	return gen_degs;
}

Poly1d Database::load_gen_diffs(const std::string& table_name) const
{
	Poly1d diffs;
	Statement stmt;
	stmt.init(*this, "SELECT gen_diff FROM " + table_name + " ORDER BY gen_id;");
	while (stmt.step() == SQLITE_ROW)
		diffs.push_back(str_to_Poly(stmt.column_str(0)));
	std::cout << "diffs loaded from " << table_name << ", size=" << diffs.size() << '\n';
	return diffs;
}

Poly1d Database::load_gen_reprs(const std::string& table_name) const
{
	Poly1d reprs;
	Statement stmt;
	stmt.init(*this, "SELECT repr FROM " + table_name + " ORDER BY gen_id;");
	while (stmt.step() == SQLITE_ROW)
		reprs.push_back(str_to_Poly(stmt.column_str(0)));
	std::cout << "reprs loaded from " << table_name << ", size=" << reprs.size() << '\n';
	return reprs;
}

Mon2d Database::load_leading_terms(const std::string& table_name) const
{
	Mon2d leadings;
	Statement stmt;
	stmt.init(*this, "SELECT leading_term FROM " + table_name + " ORDER BY t;");
	while (stmt.step() == SQLITE_ROW) {
		Mon mon(str_to_Mon(stmt.column_str(0)));
		if (size_t(mon[0].gen) >= leadings.size())
			leadings.resize(size_t(mon[0].gen) + 1);
		leadings[mon[0].gen].push_back(mon);
	}
	std::cout << "leadings loaded from " << table_name << ", size=" << leadings.size() << '\n';
	return leadings;
}

Poly1d Database::load_gb(const std::string& table_name) const
{
	Poly1d gb;
	Statement stmt;
	stmt.init(*this, "SELECT leading_term, basis FROM " + table_name + " ORDER BY t;");
	while (stmt.step() == SQLITE_ROW) {
		Mon lead(str_to_Mon(stmt.column_str(0)));
		Poly basis(str_to_Poly(stmt.column_str(1)));
		Poly g = add(basis, { lead });
		gb.push_back(std::move(g));
	}
	std::cout << "gb loaded from " << table_name << ", size=" << gb.size() << '\n';
	return gb;
}

std::map<Deg, Mon1d> Database::load_basis(const std::string& table_name) const
{
	std::map<Deg, Mon1d> basis;
	Statement stmt;
	stmt.init(*this, "SELECT s, t, v, mon FROM " + table_name + " ORDER BY mon_id;");
	while (stmt.step() == SQLITE_ROW) {
		Deg d = { stmt.column_int(0), stmt.column_int(1), stmt.column_int(2) };
		basis[d].push_back(str_to_Mon(stmt.column_str(3)));
	}
	std::cout << "basis loaded from " << table_name << ", size=" << basis.size() << '\n';
	return basis;
}

std::map<Deg, array2d> Database::load_mon_diffs_ind(const std::string& table_name) const
{
	std::map<Deg, array2d> diffs_ind;
	Statement stmt;
	stmt.init(*this, "SELECT s, t, v, diff FROM " + table_name + " WHERE diff IS NOT NULL ORDER BY mon_id;");
	while (stmt.step() == SQLITE_ROW) {
		Deg deg = { stmt.column_int(0), stmt.column_int(1), stmt.column_int(2) };
		diffs_ind[deg].push_back(str_to_array(stmt.column_str(3)));
	}
	std::cout << "diffs_ind loaded from " << table_name << ", size=" << diffs_ind.size() << '\n';
	return diffs_ind;
}

std::map<Deg, Poly1d> Database::load_mon_diffs(const std::string& table_name, const std::map<Deg, Mon1d>& basis, int r) const
{
	std::map<Deg, Poly1d> diffs;
	Statement stmt;
	stmt.init(*this, "SELECT s, t, v, diff FROM " + table_name + " WHERE diff IS NOT NULL ORDER BY mon_id;");
	while (stmt.step() == SQLITE_ROW) {
		Deg deg = { stmt.column_int(0), stmt.column_int(1), stmt.column_int(2) };
		array diff_index = str_to_array(stmt.column_str(3));
		diffs[deg].push_back(diff_index.empty() ? Poly{} : indices_to_Poly(std::move(diff_index), basis.at(deg + Deg{ 1, 0, -r })));
	}
	std::cout << "diffs loaded from " << table_name << ", size=" << diffs.size() << '\n';
	return diffs;
}

std::map<Deg, BasisComplex> Database::load_basis_ss(const std::string& table_name, int r) const
{
	std::map<Deg, BasisComplex> basis_ss;
	Statement stmt;
	stmt.init(*this, "SELECT s, t, v, level, base FROM " + table_name + " ;");
	while (stmt.step() == SQLITE_ROW) {
		Deg deg = { stmt.column_int(0), stmt.column_int(1), stmt.column_int(2) };
		int level = stmt.column_int(3);
		if (level <= r)
			basis_ss[deg].boundaries.push_back(str_to_array(stmt.column_str(4)));
		else if (level <= T_MAX - r)
			basis_ss[deg].cycles.push_back(str_to_array(stmt.column_str(4)));
	}
	std::cout << "basis_ss loaded from " << table_name << ", size=" << basis_ss.size() << '\n';
	return basis_ss;
}

void Database::save_generators(const std::string& table_name, const std::vector<Deg>& gen_degs, const Poly1d& gen_reprs, size_t i_start) const
{
	Statement stmt_update_generators;
	stmt_update_generators.init(*this, "INSERT INTO " + table_name + " (gen_id, repr, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5);");

	for (size_t i = i_start; i < gen_degs.size(); ++i) {
		stmt_update_generators.bind_int(1, (int)i);
		stmt_update_generators.bind_str(2, Poly_to_str(gen_reprs[i]));
		stmt_update_generators.bind_int(3, gen_degs[i].s);
		stmt_update_generators.bind_int(4, gen_degs[i].t);
		stmt_update_generators.bind_int(5, gen_degs[i].v);
		stmt_update_generators.step_and_reset();
	}
#ifdef DATABASE_SAVE_LOGGING
	std::cout << gen_degs.size() - i_start << " generators are inserted into " + table_name + "!\n";
#endif
}

void Database::save_gb(const std::string& table_name, const Poly1d& gb, const std::vector<Deg>& gen_degs, size_t i_start) const
{
	Statement stmt_update_relations;
	stmt_update_relations.init(*this, "INSERT INTO " + table_name + " (leading_term, basis, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5);");

	for (size_t i = i_start; i < gb.size(); ++i) {
		Deg deg = get_deg(gb[i], gen_degs);
		stmt_update_relations.bind_str(1, Mon_to_str(gb[i].front()));
		stmt_update_relations.bind_str(2, Poly_to_str(gb[i].begin() + 1, gb[i].end()));
		stmt_update_relations.bind_int(3, deg.s);
		stmt_update_relations.bind_int(4, deg.t);
		stmt_update_relations.bind_int(5, deg.v);
		stmt_update_relations.step_and_reset();
	}
#ifdef DATABASE_SAVE_LOGGING
	std::cout << gb.size() - i_start << " relations are inserted into " + table_name + "!\n";
#endif
}

void Database::save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis) const
{
	Statement stmt;
	stmt.init(*this, "INSERT INTO " + table_name + " (mon, s, t, v) VALUES (?1, ?2, ?3, ?4);");

	int count = 0;
	for (auto& [deg, basis_d] : basis) {
		for (auto& m : basis_d) {
			++count;
			stmt.bind_str(1, Mon_to_str(m));
			stmt.bind_int(2, deg.s);
			stmt.bind_int(3, deg.t);
			stmt.bind_int(4, deg.v);
			stmt.step_and_reset();
		}
	}
#ifdef DATABASE_SAVE_LOGGING
	std::cout << count << " bases are inserted into " + table_name + "!\n";
#endif
}

void Database::save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis, const std::map<Deg, array2d>& mon_reprs) const
{
	Statement stmt;
	stmt.init(*this, "INSERT INTO " + table_name + " (mon, repr, s, t, v) VALUES (?1, ?2, ?3, ?4);");

	int count = 0;
	for (auto& [deg, basis_d] : basis) {
		for (size_t i = 0; i < basis_d.size(); ++i){
			++count;
			stmt.bind_str(1, Mon_to_str(basis_d[i]));
			stmt.bind_str(2, array_to_str(mon_reprs.at(deg)[i]));
			stmt.bind_int(3, deg.s);
			stmt.bind_int(4, deg.t);
			stmt.bind_int(5, deg.v);
			stmt.step_and_reset();
		}
	}
#ifdef DATABASE_SAVE_LOGGING
	std::cout << count << " bases are inserted into " + table_name + "!\n";
#endif
}

void Database::save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis, const std::map<Deg, Poly1d>& mon_reprs) const
{
	Statement stmt;
	stmt.init(*this, "INSERT INTO " + table_name + " (mon, repr, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5);");

	for (auto& [deg, basis_d] : basis) {
		for (size_t i = 0; i < basis_d.size(); ++i) {
			stmt.bind_str(1, Mon_to_str(basis_d[i]));
			stmt.bind_str(2, Poly_to_str(mon_reprs.at(deg)[i]));
			stmt.bind_int(3, deg.s);
			stmt.bind_int(4, deg.t);
			stmt.bind_int(5, deg.v);
			stmt.step_and_reset();
		}
	}
#ifdef DATABASE_SAVE_LOGGING
	std::cout << "basis is inserted into " + table_name + ", number of degrees=" << basis.size() << '\n';
#endif
}

void Database::save_basis_ss(const std::string& table_name, const std::map<Deg, BasisSS>& basis_ss) const
{
	Statement stmt;
	stmt.init(*this, "INSERT INTO " + table_name + " (base, diff, level, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");

	for (const auto& [deg, basis_ss_d] : basis_ss) {
		for (size_t i = 0; i < basis_ss_d.basis_ind.size(); ++i) {
			stmt.bind_str(1, array_to_str(basis_ss_d.basis_ind[i]));
			stmt.bind_str(2, array_to_str(basis_ss_d.diffs_ind[i]));
			stmt.bind_int(3, basis_ss_d.levels[i]);
			stmt.bind_int(4, deg.s);
			stmt.bind_int(5, deg.t);
			stmt.bind_int(6, deg.v);
			stmt.step_and_reset();
		}
	}
#ifdef DATABASE_SAVE_LOGGING
	std::cout << "basis_ss is inserted into " + table_name + ", number of degrees=" << basis_ss.size() << '\n';
#endif
}

Statement::~Statement() { sqlite3_finalize(stmt_); }
void Statement::init(const Database& db, const std::string& sql) { db.sqlite3_prepare_v100(sql, &stmt_); }


void Statement::bind_str(int iCol, const std::string& str) const
{
	if (sqlite3_bind_text(stmt_, iCol, str.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK)
		throw "29cc3b21";
}

void Statement::bind_int(int iCol, int i) const
{
	if (sqlite3_bind_int(stmt_, iCol, i) != SQLITE_OK)
		throw "a61e05b2";
}

const char* Statement::column_str(int iCol) const
{
	return reinterpret_cast<const char*>(sqlite3_column_text(stmt_, iCol));
}

int Statement::column_int(int iCol) const
{
	return sqlite3_column_int(stmt_, iCol);
}

int Statement::column_type(int iCol) const
{
	return sqlite3_column_type(stmt_, iCol);
}

int Statement::step() const
{
	return sqlite3_step(stmt_);
}

int Statement::reset() const
{
	return sqlite3_reset(stmt_);
}

void Statement::step_and_reset() const
{
	step(); reset();
}

/* Find the index of mon in basis, assuming that mon is in basis and basis is sorted. */
inline int get_index(const Mon1d& basis, const Mon& mon)
{
	auto first = std::lower_bound(basis.begin(), basis.end(), mon);
#ifdef _DEBUG
	if (first == basis.end() || mon < (*first)) {
		std::cout << "index not found\n";
		throw "178905cf";
	}
#endif
	return int(first - basis.begin());
}

array Poly_to_indices(const Poly& poly, const Mon1d& basis)
{
	array result;
	for (const Mon& mon : poly)
		result.push_back(get_index(basis, mon));
	return result;
}

Poly indices_to_Poly(const array& indices, const Mon1d& basis)
{
	Poly result;
	for (int i : indices)
		result.push_back(basis[i]);
	return result;
}

std::string array_to_str(array::const_iterator pbegin, array::const_iterator pend)
{
	std::stringstream ss;
	for (auto p = pbegin; p < pend; p++) {
		ss << *p;
		if (p + 1 < pend)
			ss << ",";
	}
	return ss.str();
}

array str_to_array(const char* str_mon)
{
	array result;
	if (str_mon[0] == '\0')
		return array();
	std::stringstream ss(str_mon);
	while (ss.good()) {
		int i;
		ss >> i;
		result.push_back(i);
		if (ss.peek() == ',')
			ss.ignore();
	}
	return result;
}

std::string Mon_to_str(MonInd pbegin, MonInd pend)
{
	std::stringstream ss;
	for (auto p = pbegin; p < pend; p++) {
		ss << p->gen << ',' << p->exp;
		if (p + 1 < pend)
			ss << ",";
	}
	return ss.str();
}

Mon str_to_Mon(const char* str_mon)
{
	Mon result;
	if (str_mon[0] == '\0')
		return {};
	std::stringstream ss(str_mon);
	while (ss.good()) {
		int gen, exp;
		ss >> gen >> "," >> exp;
		result.emplace_back(gen, exp);
		if (ss.peek() == ',')
			ss.ignore();
	}
	return result;
}

std::string Poly_to_str(Poly::const_iterator pbegin, Poly::const_iterator pend) /* Warning: assume the algebra is connected */
{
	std::stringstream ss;
	for (auto pMon = pbegin; pMon < pend; ++pMon) {
		for (auto p = pMon->begin(); p != pMon->end(); ++p) {
			ss << p->gen << ',' << p->exp;
			if (p + 1 != pMon->end())
				ss << ",";
		}
		if (pMon + 1 != pend)
			ss << ";";
	}
	return ss.str();
}

Poly str_to_Poly(const char* str_poly) /* Warning: assume the algebra is connected */
{
	Poly result;
	if (str_poly[0] == '\0')
		return {};
	std::stringstream ss(str_poly);
	while (ss.good()) {
		int gen, exp;
		ss >> gen >> "," >> exp;
		if (result.empty())
			result.emplace_back();
		result.back().emplace_back(gen, exp);
		if (ss.peek() == ',')
			ss.ignore();
		else if (ss.peek() == ';') {
			ss.ignore();
			result.emplace_back();
		}
	}
	return result;
}